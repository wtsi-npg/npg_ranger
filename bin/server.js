#!/usr/bin/env node

"use strict";

var os      = require('os');
var fs      = require('fs');
var fse     = require('fs-extra');
var path    = require('path');
var http    = require('http');
var child     = require('child_process');
var url     = require('url');
var util    = require('util');
var MongoClient = require('mongodb').MongoClient;
var GetOpt      = require('node-getopt');

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

function createPromiseForProcess(p) {
  return new Promise(
    function(resolve, reject) {

      p.stderr.on('data', function(data) {
        console.log('STDERR FOR ' + p.title + ': ' + data);
      });

      p.on('close', function(code, signal) {
        if (code || signal) {
          reject(code || signal);
        } else {
          resolve();
        }
      });
      p.on('error', function(err) {
        reject(err);
      });
      p.stdin.on('error', (err) => {
        reject(err);
      });
      p.stdout.on('error', (err) => {
        reject(err);
      });
    }
  );
}


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
    throw 'Either some files are bam and some are cram or all files are in unexpected format';
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

  // Unfortunatelly, no way to access the trailers already set.
  // We cannot check that, if the message has gone already,
  // it had correct trailer.
  if (!response.finished) {
    var header = {};
    header[DATA_TRUNCATION_TRAILER] = success ? 'false' : 'true';
    response.addTrailers(header);
    response.end();
  }
}

function errorResponse(response, code, m) {

  m = m || 'Unknown error';
  if (code == 500) {
    m = 'Internal server error: ' + m;
  }
  console.log(m);
  if (!response.headersSent) {
    response.statusCode  = code;
    response.statusMessage = m;
  }
  response.end();
}

function getFile(response, query, user) {

  if (user) {
    console.log('IMPLEMENT USER AUTHORISATION FOR FILES FROM IRODS!');
  }

  setContentType(response, query);
  const view = child.spawn(SAMTOOLS_COMMAND, stViewAttrs(query));
  const prom = createPromiseForProcess(view);
  view.title = 'samtools view';

  prom.then(
      // Report success
      function() {
        console.log('SUCCESS ' + view.title);
        endResponse(response,1);
      }
    )
    .catch(
       // Log the rejection reason
       function(code) {
         console.log('ERROR ' + view.title + ' ' + code);
         endResponse(response,0);
       }
    );

  process.nextTick(() => {
    view.stdout.pipe(response, {end: false});
  });
}

function mergeFiles(response, query) {

  setContentType(response, query);
  var dir = tempFilePath();  console.log('DIRECTORY ' + dir);
  fs.mkdirSync(dir);

  const cleanup = function() {
    fse.remove(dir, function(err) {
      if (err) {
        console.log(`Failed to remove ${dir}: ${err}`);
      }
    });
  };

  const merge = child.spawn(SAMTOOLS_COMMAND,
                stMergeAttrs(query),
                {cwd: dir});
  const mergeprom = createPromiseForProcess(merge);
  merge.title = 'samtools merge';

  const markdup = child.spawn(BBB_MARKDUPS_COMMAND, bbbMarkDupsAttrs());
  const markdupprom = createPromiseForProcess(markdup);
  markdup.title = BBB_MARKDUPS_COMMAND;

  delete query.region;
  delete query.directory;
  const view  = child.spawn(SAMTOOLS_COMMAND, stViewAttrs(query));
  const viewprom = createPromiseForProcess(view);
  view.title = 'samtools view (post-merge)';

  var processes = [merge,markdup,view];
  var promises = [mergeprom,markdupprom,viewprom];
  promises.forEach(function(promise, i) {
    var title = processes[i].title;
    promise.then(
      // Report success
      function() {
        console.log('SUCCESS ' + title);
      }
    )
    .catch(
      // Log the rejection reason
      function(code) {
        console.log('ERROR ' + title + ' ' + code);
      }
    );
  });

  Promise.all(promises).then(
      function() { endResponse(response,1); cleanup();},
      function() { endResponse(response,0); cleanup();}
                            );

  process.nextTick(() => {
    merge.stdout.pipe(markdup.stdin);
    markdup.stdout.pipe(view.stdin);
    view.stdout.pipe(response, {end: false});
  });
}

function authorise(user, files, whatnot, badluck) {
  if (opt.options.skipauth) {
    whatnot();
  } else
  if (user.username) {
    if (files && (files instanceof Array) && files.length) {
      var agroup_ids = files.map(function(file) { return file.access_control_group_id; });
      if (agroup_ids.every(function(id) { return id; })) {
        agroup_ids.sort();
        // Get a list of unique ids
        agroup_ids = agroup_ids.filter(function(item, index, thisArray) {
          return (index === 0) ? 1 : ((item === thisArray[index - 1]) ? 0 : 1);
        });
        console.log('ACCESS GROUP IDS: ' + agroup_ids.join(' '));
        var dbquery = {
          members:                 user.username,
          access_control_group_id: {$in: agroup_ids}
        };
        db.collection('access_control_group').count(dbquery,
          function(err, count) {
            if (err) {
              badluck('Failed to get authorisation info');
              return;
            }
            if (count == agroup_ids.length) {
              whatnot();
            } else {
              var qualifier = count ? 'any' : 'some';
              badluck('Not authorised for ' + qualifier + ' of the files');
            }
          }
        );
      } else {
        badluck('Access group id is missing for one of the files');
      }
    } else {
      badluck('File info is not available, cannot authorise access');
    }
  } else {
    badluck('Username is not known');
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
    throw 'Sample accession number or file should be given';
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
      query.files = files.map(function(f) { return f.filepath_by_host[HOST] || f.filepath_by_host["*"]; });
      var numFiles = files.length;
      if (numFiles === 0) {
        console.log('No files for ' + (a ? a : query.name) );
        response.end();
      } else {
        var whatnot = function() {
          console.log("User " + user.username + " is given access");
          if (numFiles === 1) {
            getFile(response, query);
          } else {
            mergeFiles(response, query);
          }
        };
        var badluck = function(message) {
          var m = util.format(
            "Authorisation failed for %s: %s", user.username, message);
          console.log(m);
          errorResponse(response, 401, m);
        };

        authorise(user, files, whatnot, badluck);
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

  try {
    var url_obj = url.parse(request.url, true);
    var path = url_obj.pathname;
    path = path ? path : '';
    var q = url_obj.query;

    if (!q.format) {
      q.format = DEFAULT_FORMAT;
    }

    response.setHeader('Trailer', DATA_TRUNCATION_TRAILER);
    var user = getUser(request);

    switch (path) {
      case '/file': {
        getData(response, q, user);
        break;
      }
      case '/sample': {
        getData(response, q, user);
        break;
      }
      default: {
        errorResponse(response, 404, 'Not found: ' + request.url);
      }
    }

  } catch (err) {
    errorResponse(response, 500,
      'Error handling request for ' + request.url + ': ' + err);
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
