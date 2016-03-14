#!/usr/bin/env node

var os          = require('os');
var fs          = require('fs');
var path        = require('path');
var http        = require('http');
var child       = require('child_process');
var url         = require('url');
var MongoClient = require('mongodb').MongoClient;
var crypto = require('crypto');

const MONGO                = 'mongodb://sf2-farm-srv1:27017/imetacache';
const SAMTOOLS_COMMAND     = 'samtools';
const BBB_MARKDUPS_COMMAND = 'bamstreamingmarkduplicates';
const IRODS_PATH_PREFIX    = 'irods:';
const TEMP_DATA_DIR_NAME   = 'npg_ranger_data';
const TEMP_DATA_DIR        = path.join(os.tmpdir(), process.env.USER, TEMP_DATA_DIR_NAME);

var db;

function socket_path() {
    var howMany = 10;
    var chars = "abcdefghijklmnopqrstuwxyzABCDEFGHIJKLMNOPQRSTUWXYZ0123456789";
    var len   = chars.length;
    var value = value = new Array(howMany);
    var rnd   = crypto.randomBytes(howMany);

    for (var i = 0; i < howMany; i++) {
        value[i] = chars[rnd[i] % len]
    };

    return '/tmp/' + value.join('') + '/sock';
}

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
    
    var file = '-';
    if (query.directory && query.name) {
        file = query.directory + '/' + query.name;
        if (query.irods) {
            file = IRODS_PATH_PREFIX + file;
        }
    }
    attrs.push(file);

    if (query.region) {
        attrs = attrs.concat(query.region);
    }
    
    console.log(attrs);
    return attrs;
}

function stMergeAttrs(query) {

    var attrs = ['merge'];
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

    attrs = attrs.concat(files.map(function(f){
        return IRODS_PATH_PREFIX + f.collection + '/' + f.data_object;
    }));

    console.log(attrs);
    return attrs;
}

function bbbMarkDupsAttrs() {
    var attrs = ['level=0','verbose=0','resetdupflag=1'];
    attrs.push('tmpfile=' + tempFilePath());
    attrs.push('M=' + tempFilePath());
    return attrs;
}

function setProcessCallbacks (pr, child, response) {

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
        var m; 
        if (code) {
            m = title + ' exited (on exit) with code ' + code;
            console.log(m);
            if (child) {
                child.kill();
            } else {
                errorResponse(response, 500, m);
            }
        }
    });

    pr.on('close', function (code, signal) {
        if (code) {
            var m = title + ' exited (on close) with code ' + code;
            if (!child) {
                errorResponse(response, 500, m);
	    }
        } else if (signal != null) {
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
    } else {
        response.addTrailers({'data_truncated': 'true'});
    }
    response.end();
}

function getFile(response, query){

    if (!(query.directory && query.name)) {
        throw 'Both directory and name should be given';
    }
    setContentType(response, query);
    const view = child.spawn(SAMTOOLS_COMMAND, stViewAttrs(query));
    view.title = 'samtools view';
    view.stdout.pipe(response);
    setProcessCallbacks(view, null, response);
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

    setProcessCallbacks(merge,   markdup, response);
    setProcessCallbacks(markdup, view,    response);
    setProcessCallbacks(view,    null,    response);
}

function getSampleData(response, query){

    var a = query.accession;
    if (!a) {
        throw 'Sample accession number should be given';
    }
    
    var dbquery = { $and: [{'avh.sample_accession_number': a},
                           {'avh.target':    "1"},
                           {'avh.manual_qc': "1"},
                           {'avh.alignment': "1"} ]};
    var columns = {_id:0, collection:1, data_object: 1};
    var files   = [];
    
    var cursor = db.collection('fileinfo').find(dbquery, columns);
    cursor.each(function(err, doc) {
        if(err) throw err;
        if (doc != null) {
            files.push(doc);
        } else {
            var numFiles = files.length;
            if (numFiles === 0) {
                console.log('No files for sample accession ' + a);
                response.end();
            } else if (numFiles === 1) {
                var d = files[0];
                query.directory = d.collection;
                query.name      = d.data_object;
                query.irods     = 1;
                getFile(response, query);
            } else {
                query.files = files;
                mergeFiles(response, query);
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

function handleRequest(request, response){

    try {
        var url_obj = url.parse(request.url, true);
        var path = url_obj.pathname;
        path = path ? path : '';
        var q = url_obj.query;

        if (!q.format) {
            q.format = 'bam';
        }

        response.setHeader('Trailer', 'data_truncated'); 

    	switch(path) {
            case '/test':
            runTest(request, response, q);
                break;
            case '/file':
                getFile(response, q);
                break;
            case '/sample':
                getSampleData(response, q);
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
    console.log('Connected to mongodb');

    //Callback for a graceful exit
    server.on('close', function(){
         console.log('Database connection closing');
         database.close();
         console.log('Server closing');
    });

    createTempDataDir();
    
    var sock = process.argv[2] || '/tmp/' + process.env.USER + '/npg_ranger.sock';
    //Lets start our server
    server.listen(sock, function(){
        //Callback triggered when server is successfully listening. Hurray!
        console.log("Server listening on %s, %s", os.hostname(), sock);
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
