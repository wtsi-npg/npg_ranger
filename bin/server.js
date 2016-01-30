var http  = require('http');
var child = require('child_process');
var url   = require('url');
var MongoClient = require('mongodb').MongoClient;

const PORT=9444;
const MONGO='mongodb://sf-nfs-01-01:27017/imetacache';

var db;

function getFile(response, query){

    var file = query.directory + '/' + query.name;
    if (query.irods) {
        file = 'irods:' + file;
    }
    var attrs = ['view', '-h'];
    if (query.format && (query.format === 'bam' || query.format === 'cram')) {
        response.setHeader("Content-Type", 'application/octet-stream');
        attrs.push(query.format === 'bam' ? '-b' : '-C');
    }
    attrs.push(file);
    if (query.region) {
        attrs = attrs.concat(query.region);
    }
    
    console.log(attrs);
    const bam = child.spawn('samtools1', attrs);
    bam.stdout.pipe(response);
    
    bam.stdout.on('end', function () {
        console.log('body finished');
    });
    bam.stderr.on('data', function (data) {
        console.log(data);
    });

    bam.on('close', function (code) {
        // Would be good to have the error itself
        console.log('child process exited with code ' + code);
    });
}

function mergeSample(response, query){

    var a = query.accession;
    if (!a) {
        throw 'Sample accession number should be given';
    }
    
    var query   = { $and: [
        {avus:{$elemMatch:{attribute: 'sample_accession_number',value: a}}},
        {avus:{$elemMatch:{attribute: 'target'                 ,value: "1"}}},
        {avus:{$elemMatch:{attribute: 'manual_qc'              ,value: "1"}}},
        {avus:{$elemMatch:{attribute: 'alignment'              ,value: "1"}}} ]};
    var columns = {_id:0, collection:1, data_object: 1};
    var cursor = db.collection('fileinfo').find(query, columns);
    var files = [];
    cursor.each(function(err, doc) {
        if(err) throw err;
        if (doc != null) {
            files.push(doc.collection + '/' + doc.data_object);
        } else {
            //file merge should go here
            response.end("Files to merge: " + files + "\n");
        }
    });
}

function handleRequest(request, response){

    try {
        var url_obj = url.parse(request.url, true);
        var path = url_obj.pathname;
        path = path ? path : '';
        var q = url_obj.query;

    	switch(path) {
            case '/file':
                getFile(response, q);
                break;
            case '/sample':
                mergeSample(response, q);
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
    server.listen(PORT, function(){
        //Callback triggered when server is successfully listening. Hurray!
        console.log("Server listening on: http://localhost:%s", PORT);
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
