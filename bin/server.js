

//Lets require/import the HTTP module
var http = require('http');
var child = require('child_process');
var url = require('url');

//x  var   fs = require('fs');

//Lets define a port we want to listen to
const PORT=9444;

function get_file(response, query){

    console.log(query.name + ' ' + query.directory + ' ' +  query.region);
    var file = query.directory + '/' + query.name;
    const bam_header = child.spawn('samtools', ['view', '-H', file] );
    //const bam        = child.spawn('samtools', ['view', file, query.region] );
    const bam        = child.spawn('samtools', ['view', file] );
    bam_header.stdout.pipe(response, { end: false });
    bam_header.stdout.on('end', function () {
       console.log('header finished');
	//return false;
    })
    
    bam_header.on('close', function (code) {
      console.log('child process exited with code ' + code);
      //return false;
    });
    
  
    //response.setHeader("Content-Type", 'application/octet-stream');
    var error = '';
    bam.stdout.pipe(response);
    
    //bam.stdout.on('data', function (data) {
      
      //console.log('stdout: ' + data);
      //console.log('piping');
      //bam.stdout.pipe(response);
    //});
    bam.stdout.on('end', function () {
       console.log('body finished');
    });
     bam.stderr.on('data', function (data) {
       console.log(data);
    });

    bam.on('close', function (code) {
      console.log('child process exited with code ' + code);
    });
}

//We need a function which handles requests and send response
function handleRequest(request, response){
    try {
        var url_obj = url.parse(request.url, true);
        var path = url_obj.pathname;
        path = path ? path : '';
        console.log(path);
        var q = url_obj.query;

	switch(path) {
            case '/file':
                console.log('case file');
                get_file(response, q)
                break;
            default:
                response.statusCode = 404;
                response.end('Not found: ' + request.url);
	}

    //response.end('It Works!! Path Hit: ' + request.url);
    //request.pipe(bam);
    } catch (err) {
        console.log(err);
        response.statusCode = 500;
        response.end('Error: ' + err);
    }
}

//Create a server
var server = http.createServer(handleRequest);

//Lets start our server
server.listen(PORT, function(){
    //Callback triggered when server is successfully listening. Hurray!
    console.log("Server listening on: http://localhost:%s", PORT);
});