

//Lets require/import the HTTP module
var http = require('http');
var child = require('child_process');
var url = require('url');

//x  var   fs = require('fs');

//Lets define a port we want to listen to
const PORT=9444;

function get_file(response, query){

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

//We need a function which handles requests and send response
function handleRequest(request, response){
    console.log('handling request');
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

//Lets start our server
server.listen(PORT, function(){
    //Callback triggered when server is successfully listening. Hurray!
    console.log("Server listening on: http://localhost:%s", PORT);
});