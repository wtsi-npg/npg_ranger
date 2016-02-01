#!/usr/bin/env node

var http        = require('http');
var child       = require('child_process');
var url         = require('url');
var MongoClient = require('mongodb').MongoClient;

const PORT                = 9444;
const MONGO               = 'mongodb://sf2-farm-srv1:27017/imetacache';
const SAMTOOLS_COMMAND    = 'samtools';
const IRODS_PATH_PREFIX   = 'irods:';

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

function getFile(response, query){

    setContentType(response, query);

    const view = child.spawn(SAMTOOLS_COMMAND, stViewAttrs(query));
    view.stdout.pipe(response);
    
    view.stdout.on('end', function () {
        console.log('Samtools view finished');
    });
    view.on('exit', function (code) {
        var err = view.stderr.read();
        if (err) {
            console.log('Samtools view stderr: ' + err);
        }
    });
    view.on('close', function (code) {
        console.log('Samtools view exited with code ' + code);
    });
}

function mergeFiles(response, query){

    setContentType(response, query);

    const merge = child.spawn(SAMTOOLS_COMMAND, stMergeAttrs(query));
    delete query['region'];
    delete query['directory'];
    //const view  = child.spawn(SAMTOOLS_COMMAND, stViewAttrs(query));
    //merge.stdout.pipe(view);
    //view.stdout.pipe(response);
    merge.stdout.pipe(response);

    merge.stdout.on('end', function () {
        console.log('Samtools merge finished');
    });
    merge.on('exit', function (code) {
        if (code) {
            console.log('Child process stderr: ' + (merge.stderr.read() || ''));
        }
    });
    merge.on('close', function (code) {
        console.log('Child process exited with code ' + code);
    });

    /* view.stdout.on('end', function () {
        console.log('Samtools view finished');
    });
    view.on('exit', function (code) {
        var err = view.stderr.read();
        if (err) {
            console.log('Child process stderr: ' + err);
        }
    });
    view.on('close', function (code) {
        console.log('Child process exited with code ' + code);
    }); */
}

function getSampleData(response, query){

    var a = query.accession;
    if (!a) {
        throw 'Sample accession number should be given';
    }
    
    var dbquery   = { $and: [
        {avus:{$elemMatch:{attribute: 'sample_accession_number',value: a}}},
        {avus:{$elemMatch:{attribute: 'target'                 ,value: "1"}}},
        {avus:{$elemMatch:{attribute: 'manual_qc'              ,value: "1"}}},
        {avus:{$elemMatch:{attribute: 'alignment'              ,value: "1"}}} ]};
    var columns = {_id:0, collection:1, data_object: 1};
    var files = [];
    
    var cursor = db.collection('fileinfo').find(dbquery, columns);
    cursor.each(function(err, doc) {
        if(err) throw err;
        if (doc != null) {
            files.push(doc);
        } else {
            var numFiles = files.length;
            if (numFiles == 0) {
                console.log('No files for sample accession ' + a);
                response.end();
            } else if (numFiles == 1) {
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

function handleRequest(request, response){

    try {
        var url_obj = url.parse(request.url, true);
        var path = url_obj.pathname;
        path = path ? path : '';
        var q = url_obj.query;

        if (!q.format) {
            q.format = 'bam';
        }

    	switch(path) {
            case '/file':
                getFile(response, q);
                break;
            case '/sample':
                getSampleData(response, q);
                break;
            default:
                response.statusCode = 404;
                var m = 'Not found: ' + request.url;
                console.log(m);
                response.end(m);
    	}

    } catch (err) {
        var m = 'Error: ' + err;
        console.log(m);
        response.statusCode = 500;
        response.end(m);
    }
}

var customPort = process.argv[2] || PORT;

//Create a server
var server = http.createServer(handleRequest);

var mongo_options = {
  db:{
    numberOfRetries : 5
  },
  server: {
    auto_reconnect: true,
    poolSize : 40,
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

    //Lets start our server
    server.listen(customPort, function(){
        //Callback triggered when server is successfully listening. Hurray!
        console.log("Server listening on: http://localhost:%s", customPort);
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
