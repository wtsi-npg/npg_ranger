#!/usr/bin/env node

var os          = require('os');
var fs          = require('fs');
var path        = require('path');
var http        = require('http');
var child       = require('child_process');
var url         = require('url');
var util        = require('util');
var MongoClient = require('mongodb').MongoClient;
var GetOpt      = require('node-getopt');

var opt = new GetOpt([
    ['p','port=PORT'        ,'PORT or socket server listens on'],
    ['m','mongourl=URI'     ,'URI to contact mongodb'],
    ['t','tempdir=PATH'     ,'PATH of temporary directory'],
    ['H','hostname=HOST'    ,'override hostname with HOST'],
    ['s','skipauth'         ,'skip authorisation steps'],
    ['h','help'             ,'display this help']
]).bindHelp().parseSystem();

const PORT                 = opt.options.port || opt.argv[0] || 9444;
const HOST                 = opt.options.hostname || os.hostname() || 'localhost';
const MONGO                = opt.options.mongourl || 'mongodb://sf2-farm-srv1:27017/imetacache';
const SAMTOOLS_COMMAND     = 'samtools';
const BBB_MARKDUPS_COMMAND = 'bamstreamingmarkduplicates';
const TEMP_DATA_DIR_NAME   = 'npg_ranger_data';
const DATA_TRUNCATION_TRAILER = 'data-truncated';
const TEMP_DATA_DIR        = opt.options.tempdir || path.join(os.tmpdir(), process.env.USER, TEMP_DATA_DIR_NAME);

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
    
    console.log(attrs);
    return attrs;
}

function stMergeAttrs(query) {

    var attrs = ['merge', "-u"];
    if (query.region) {
        var regions = query.region;
        if (typeof regions != 'object') {
            regions = [regions];
        }
        regions.map(function(r){
            attrs.push('-R');
            attrs.push(r);
        });
    }
    attrs.push('-');
 
    var files = query.files;
    var re_bam  = /\.bam$/;
    var re_cram = /\.cram$/;
    var some_bam  = files.some(function(f){ return  re_bam.test(f.data_object)});
    var some_cram = files.some(function(f){ return re_cram.test(f.data_object)});
    if (some_bam && some_cram) {
        throw 'Either some files are bam and some are cram or all files are in unexpected format';
    }

    attrs = attrs.concat(files );
    query.files.length=0;

    console.log(attrs);
    return attrs;
}

function bbbMarkDupsAttrs() {
    var attrs = ['level=0','verbose=0','resetdupflag=1'];
    attrs.push('tmpfile=' + tempFilePath());
    attrs.push('M=/dev/null');
    return attrs;
}

function setProcessCallbacks (pr, isHead, child, response) {

    var title = pr.title;   

    pr.stderr.on('data', function (data) {
        console.log('STDERR for ' + title + ': ' + data);
    });

    pr.on('error', function (err) {
        var m =  'Error creating process ' + title + ' ' + err;
        if (child) {
            child.kill();
        }
        errorResponse(response, 500, m);
    });

    pr.on('exit', function (code) {
        if (code) {
            if (child) {
                child.kill();
            } 
            errorResponse(response, 500,
                title + ' exited (on exit) with code ' + code);
        }
    });

    pr.on('close', function (code, signal) {
        if (code) {  
            errorResponse(response, 500,
                title + ' exited (on close) with code ' + code);  
        }  else if (signal != null) {
            console.log(title + ' terminated by a parent ' + signal);
            if (child) {
                child.kill();
	    }
	}
    });
}

function errorResponse (response, code, m) {

    m = m || '';
    if (code == 500) {
        m = m ? ('Internal server error: ' + m) : 'Internal server error';
    }
    console.log(m);
    if (!response.headersSent) {
        response.statusCode    = code;
        response.statusMessage = m;
    }
    response.end();
}

function getFile(response, query, authorised){
  
    if (!authorised) {
        throw 'Authorisation flag is not set';
    }
    setContentType(response, query);
    const view = child.spawn(SAMTOOLS_COMMAND, stViewAttrs(query));
    view.title = 'samtools view';
    view.stdout.pipe(response);
    setProcessCallbacks(view, 1, null, response);
}

function mergeFiles(response, query){

    setContentType(response, query);

    const merge = child.spawn(SAMTOOLS_COMMAND, stMergeAttrs(query));
    merge.title = 'samtools merge';

    const markdup = child.spawn(BBB_MARKDUPS_COMMAND, bbbMarkDupsAttrs());
    markdup.title = BBB_MARKDUPS_COMMAND;

    delete query['region'];
    delete query['directory'];
    const view  = child.spawn(SAMTOOLS_COMMAND, stViewAttrs(query));
    view.title = 'samtools view (post-merge)';

    merge.stdout.pipe(markdup.stdin);
    markdup.stdout.pipe(view.stdin);
    view.stdout.pipe(response);

    setProcessCallbacks(merge,   1, markdup, response);
    setProcessCallbacks(markdup, 0, view,    response);
    setProcessCallbacks(view,    0, null,    response);
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
                    return (index == 0) ? 1 : ((item === thisArray[index-1]) ? 0 : 1)});
                console.log('ACCESS GROUP IDS: ' + agroup_ids.join(' '));
                var dbquery = {
                    "members"                : user.username,
                    "access_control_group_id": {$in: agroup_ids}
                              };
                db.collection('access_control_group').count(dbquery,
                    function (err, count) {
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

function getSampleData(response, query, user){

    var a = query.accession;
    if (!a) {
        throw 'Sample accession number should be given';
    }
    
    var dbquery = { $and: [{'avh.sample_accession_number': a},
                           {'avh.target':    "1"},
                           {'avh.manual_qc': "1"},
                           {'avh.alignment': "1"} ]};
    var localkey = 'filepath_by_host.' + HOST;
    var columns = {_id:0, 'filepath_by_host.*':1, 'access_control_group_id': 1};
    columns[localkey]=1;
    var files   = [];
    
    var cursor = db.collection('fileinfo').find(dbquery, columns);
    cursor.each(function(err, doc) {
        if(err) throw err;
        if (doc != null) {
            files.push(doc);
        } else {
            cursor.close(); // Got all results, do not need the cursor any longer.
            query.files = files.map(function(f){ return f.filepath_by_host[HOST] || f.filepath_by_host["*"]; });
            var numFiles = files.length;
            if (numFiles === 0) {
                console.log('No files for sample accession ' + a);
                response.end();
            } else {
                var whatnot = function() {
                    console.log("User " + user.username + " is given access");
                    if (numFiles === 1) {
                        getFile(response, query, 1);
                    } else {
                        mergeFiles(response, query);
                    }
		}
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

function runTest(request, response, q) {

    console.log('====RAW HEADERS ' + request.rawHeaders);
    // The headers object contains 'normalized' headers,
    // note change of case.
    var user = request.headers['x-remote-user'] || 'Unknown';
    console.log('====REMOTE USER IS ' + user);

    response.setHeader("Content-Type", 'text/html');
    response.write('<html><body><h1>Hello, ' +  user + '</h1></body></html>');
    response.end();
}

function getUser(request) {
    var user = {};
    user.username   = request.headers['x-remote-user'] || null;
    return user;
}

function handleRequest(request, response){

    try {
        var url_obj = url.parse(request.url, true);
        var path = url_obj.pathname;
        path = path ? path : '';
        var q = url_obj.query;

        if (!q.format) {
            q.format = 'bam';
        }

        response.setHeader('Trailer', DATA_TRUNCATION_TRAILER); 

    	switch(path) {
            case '/test':
                runTest(request, response, q);
                break;
            case '/file':
                getFile(response, q);
                break;
            case '/sample':
                getSampleData(response, q, getUser(request));
                break;
            default:
                errorResponse(response, 404, 'Not found: ' + request.url);
    	}

    } catch (err) {
        console.log('Error handling request for ' + request.url + ': ' + err);
        errorResponse(response, 500);
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

//Create a server
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

    if(err) throw err;
    db = database;
    console.log('Connected to mongodb at ' + MONGO);

    //Callback for a graceful exit
    server.on('close', function(){
         console.log('Database connection closing');
         database.close();
         console.log('Server closing');
    });

    createTempDataDir();
    
    var sock = process.argv[2] || '/tmp/' + process.env.USER + '/npg_ranger.sock';
    //Lets start our server
    server.listen(PORT, function(){
        //Callback triggered when server is successfully listening. Hurray!
        console.log("Server listening on: http://%s:%s", HOST, PORT);
    });
});

process.on('SIGTERM', function () {
    server.close(function () {
        process.exit(0);
    });
});

process.on('SIGINT', function () {
    server.close(function () {
        process.exit(0);
    });
});
