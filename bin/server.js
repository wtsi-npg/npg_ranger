#!/usr/bin/env node

"use strict";

var os      = require('os');
var fs      = require('fs');
var fse     = require('fs-extra');
var path    = require('path');
var http    = require('http');
var child   = require('child_process');
var url     = require('url');
var MongoClient = require('mongodb').MongoClient;
var GetOpt      = require('node-getopt');

var pipeline    = require('../lib/pipeline.js');
var DataAccess  = require('../lib/auth.js');

var opt = new GetOpt([
    ['p','port=PORT'        ,'PORT or socket which server listens on'],
    ['m','mongourl=URI'     ,'URI to contact mongodb'],
    ['t','tempdir=PATH'     ,'PATH of temporary directory'],
    ['H','hostname=HOST'    ,'override hostname with HOST'],
    ['s','skipauth'         ,'skip authorisation steps'],
    ['h','help'             ,'display this help']
]).bindHelp().parseSystem();

const PORT                    = opt.options.port || opt.argv[0]
                                || path.join(os.tmpdir(), process.env.USER, 'npg_ranger.sock');
const HOST                    = opt.options.hostname || os.hostname() || 'localhost';
const MONGO                   = opt.options.mongourl || 'mongodb://sf2-farm-srv1:27017/imetacache';
const SAMTOOLS_COMMAND        = 'samtools';
const BBB_MARKDUPS_COMMAND    = 'bamstreamingmarkduplicates';
const TEMP_DATA_DIR_NAME      = 'npg_ranger_data';
const TEMP_DATA_DIR           = opt.options.tempdir || path.join(os.tmpdir(), process.env.USER, TEMP_DATA_DIR_NAME);
const DATA_TRUNCATION_TRAILER = 'data-truncated';
const DEFAULT_FORMAT          = 'bam';

var db;

function setContentType(response, query) {
  if (query.format && (query.format === 'bam' || query.format === 'cram')) {
    response.setHeader("Content-Type", 'application/octet-stream');
  }
}

function stViewAttrs(query) {
  var attrs = ['view', '-h'];

  if (query.format && (query.format === 'bam' || query.format === 'cram')) {
    attrs.push(query.format === 'bam' ? '-b' : '-C');
  }
  attrs.push(query.files.shift() || "-");
  if (query.region) {
    attrs = attrs.concat(query.region);
  }
  console.log('view attrs: ' + attrs);
  return attrs;
}

function stMergeAttrs(query) {
  var attrs = ['merge', "-u"];
  if (query.region) {
    var regions = query.region;
    if (typeof regions != 'object') {
      regions = [regions];
    }
    regions.map(function(r) {
      attrs.push('-R');
      attrs.push(r);
    });
  }
  attrs.push('-');

  var files = query.files;
  var re_bam  = /\.bam$/;
  var re_cram = /\.cram$/;
  var some_bam  = files.some(function(f) { return re_bam.test(f.data_object); });
  var some_cram = files.some(function(f) { return re_cram.test(f.data_object); });
  if (some_bam && some_cram) {
    throw new Error(
      'Inconsistent format, all files should be either bam or cram');
  }

  attrs = attrs.concat(files );
  query.files.length = 0;

  console.log('merge attrs: ' + attrs);
  return attrs;
}

function bbbMarkDupsAttrs() {
  var attrs = ['level=0','verbose=0','resetdupflag=1'];
  attrs.push('tmpfile=' + tempFilePath());
  attrs.push('M=/dev/null');
  return attrs;
}

function endResponse(response, success) {
  if (!response.finished) {
    var header = {};
    header[DATA_TRUNCATION_TRAILER] = success ? 'false' : 'true';
    response.addTrailers(header);
    response.end();
  }
}

function errorResponse(response, code, m) {
  if (!code) {
    throw new ReferenceError('Error code required');
  }
  m = m || 'Unknown error';
  if (!response.headersSent) {
    response.statusCode    = code;
    response.statusMessage = m;
  }
  response.end();
  console.log(`Server error ${code}: ${m}.`);
}

function getFile(response, query) {
  const view = child.spawn(SAMTOOLS_COMMAND, stViewAttrs(query));
  view.title = 'samtools view';
  pipeline(
    [view],
    () => {endResponse(response,1);},
    () => {endResponse(response,0);} )
    .run(response);
}

function mergeFiles(response, query) {

  var dir = tempFilePath();
  fs.mkdirSync(dir);

  const cleanup = function() {
    fse.remove(dir, function(err) {
      if (err) {
        console.log(`Failed to remove ${dir}: ${err}`);
      }
    });
  };

  var attrs;
  try {
    attrs = stMergeAttrs(query);
  } catch (ex) {
    errorResponse(response, 500, ex);
  }
  if (!attrs) {
    return;
  }

  const merge = child.spawn(SAMTOOLS_COMMAND,
                attrs,
                {cwd: dir});
  merge.title = 'samtools merge';

  const markdup = child.spawn(BBB_MARKDUPS_COMMAND, bbbMarkDupsAttrs());
  markdup.title = BBB_MARKDUPS_COMMAND;

  delete query.region;
  delete query.directory;
  const view  = child.spawn(SAMTOOLS_COMMAND, stViewAttrs(query));
  view.title = 'samtools view (post-merge)';

  pipeline(
    [merge,markdup,view],
    () => {endResponse(response,1); cleanup();},
    () => {endResponse(response,0); cleanup();} )
    .run(response);
}

function setupPipeline(response, query) {
  if (!query.format) {
    query.format = DEFAULT_FORMAT;
  }
  response.setHeader('Trailer', DATA_TRUNCATION_TRAILER);
  setContentType(response, query);
  if (query.files.length === 1) {
    getFile(response, query);
  } else {
    mergeFiles(response, query);
  }
}

function getData(response, query, user) {

  var a = query.accession;
  var dbquery;
  if (a) {
    dbquery =  { $and: [{'avh.sample_accession_number': a},
               {'avh.target':    "1"},
               {'avh.manual_qc': "1"},
               {'avh.alignment': "1"} ]};
  } else if (query.name) {
    dbquery =  {data_object: query.name};
    if (query.directory) {
      dbquery =  { $and: [ dbquery, {collection: query.directory} ]};
    }
  } else {
    throw new Error('Sample accession number or file should be given');
  }
  console.log(dbquery);

  var localkey = 'filepath_by_host.' + HOST;
  var columns = {_id: 0, 'filepath_by_host.*': 1, access_control_group_id: 1};
  columns[localkey] = 1;
  var files   = [];

  var cursor = db.collection('fileinfo').find(dbquery, columns);
  cursor.each(function(err, doc) {
    if (err) {
      throw err;
    }
    if (doc != null) {
      files.push(doc);
    } else {
      cursor.close(); // Got all results, do not need the cursor any longer.
      if (files.length === 0) {
        errorResponse(response, 404, 'No files for ' + (a ? a : query.name));
      } else {
        query.files = files.map(function(f) { return f.filepath_by_host[HOST] || f.filepath_by_host["*"]; });
        if (opt.options.skipauth) {
          setupPipeline(response, query);
        } else {
          var da = new DataAccess(db, files);
          da.on('authorised', (username) => {
            console.log(`User ${username} is given access`);
            setupPipeline(response, query);
          });
          da.on('failed', (username, message) => {
            errorResponse(response, 401,
              `Authorisation failed for user '${username}': ${message}`);
          });
          da.authorise(user.username);
        }
      }
    }
  });
}

function getUser(request) {
  var user = {};
  user.username   = request.headers['x-remote-user'] || null;
  return user;
}

function handleRequest(request, response) {

  var user = getUser(request);
  if (!user.username && !opt.options.skipauth) {
    return errorResponse(response, 407, 'Proxy authentication required');
  }

  var url_obj = url.parse(request.url, true);
  var path = url_obj.pathname;
  path = path ? path : '';
  var q = url_obj.query;

  switch (path) {
    case '/file': {
      if (!q.name) {
        errorResponse(response, 400,
          'Invalid request: file name should be given');
      } else {
        getData(response, q, user);
      }
      break;
    }
    case '/sample': {
      if (!q.accession) {
        errorResponse(response, 400,
          'Invalid request: sample accession number should be given');
      } else {
        getData(response, q, user);
      }
      break;
    }
    default: {
      errorResponse(response, 400, 'URL not available: ' + path);
    }
  }
}

function tempFilePath() {
  return path.join(TEMP_DATA_DIR, Math.random().toString().substr(2));
}

function createTempDataDir() {
  if (!fs.existsSync(TEMP_DATA_DIR)) {
    var dir = path.join(os.tmpdir(), process.env.USER);
    if (!fs.existsSync(dir)) {
      fs.mkdirSync(dir);
    }
    fs.mkdirSync(TEMP_DATA_DIR);
    console.log('Created temp data directory ' + TEMP_DATA_DIR);
  } else {
    console.log('Found temp data directory ' + TEMP_DATA_DIR);
  }
}

// Create a server
var server = http.createServer(handleRequest);

var mongo_options = {
  db: {
    numberOfRetries: 5
  },
  server: {
    auto_reconnect: true,
    poolSize: 40,
    socketOptions: {
      connectTimeoutMS: 5000
    }
  },
  replSet: {},
  mongos: {}
};

MongoClient.connect(MONGO, mongo_options, function(err, database) {

  if (err) {
    throw err;
  }
  db = database;
  console.log('Connected to mongodb at ' + MONGO);

  // Callback for a graceful exit
  server.on('close', function() {
    console.log('Database connection closing');
    database.close();
    console.log('Server closing');
  });

  createTempDataDir();

  // Lets start our server
  server.listen(PORT, function() {
    // Callback triggered when server is successfully listening. Hurray!
    console.log("Server listening on %s, %s", HOST, PORT);
  });
});

process.on('SIGTERM', function() {
  server.close(function() {
    process.exit(0);
  });
});

process.on('SIGINT', function() {
  server.close(function() {
    process.exit(0);
  });
});
