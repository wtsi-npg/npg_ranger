#!/usr/bin/env node
var http        = require('http');

var user = process.env.USER;
http.get({
       socketPath: '/tmp/' + user + '/npg_ranger.sock',
       path: '/test',
       headers: {'X-Remote-User': user}
    }, function(response) {
        // Continuously update stream with data
        var body = '';
        response.on('data', function(d) {
            body += d;
        });
        response.on('end', function() {
            // Data reception is done, do whatever with it!
            console.log("RECEIVED: " + body);
        });
 });