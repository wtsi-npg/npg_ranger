

//Lets require/import the HTTP module
var http  = require('http');
var child = require('child_process');
var url   = require('url');
var MongoClient = require('mongodb').MongoClient;

//Lets define a port we want to listen to
const PORT=9444;
const MONGO='mongodb://sf-nfs-01-01:27017/imetacache';

var db;

function getFile(response, query){

    //console.log(query.name + ' ' + query.directory + ' ' +  query.region);
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

function mergeSample(response, query, db){

    if (!query.name && !query.accession) {
        throw 'Either sample name or accession should be given';
    }
    response.end('Merging files by sample accession - wait for implementation');    
}

function handleRequest(request, response){

    db.collection("fileinfo").find({$and:[{"avus.attribute" : "library"},{"avus.value" : "11144796"}]}, function(err, docs) {
        if(err) throw err;
        console.log('some reply');
	docs.each(function(err, doc) {
            if(err) throw err;
            if(doc) {
                console.log(doc);
            } else {console.log('nothing');}
        });
    });
    try {
        var url_obj = url.parse(request.url, true);
        var path = url_obj.pathname;
        path = path ? path : '';
        console.log(path);
        var q = url_obj.query;

    	switch(path) {
            case '/file':
                getFile(response, q);
                break;
            case '/sample':
                mergeSample(response, q, db);
                break;
            default:
                response.statusCode = 404;
                console.log('Not found: ' + request.url);
                response.end('Not found: ' + request.url);
    	}

    } catch (err) {
        console.log(err);
        response.statusCode = 500;
        response.end('Error: ' + err);
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
  console.log('Connected to mongo');

  //Lets start our server
  server.listen(PORT, function(){
    //Callback triggered when server is successfully listening. Hurray!
    console.log("Server listening on: http://localhost:%s", PORT);
  });
});