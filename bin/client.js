#!/usr/bin/env node

"use strict";

var http    = require('http');

// HTTP client to connect to an HTTP server runnign on a
// local unix socket.

var user = process.env.USER;
http.get(
  {
    socketPath: '/tmp/' + user + '/npg_ranger.sock',
    path: '/sample?region=1:77970-77980&accession=ERS1060068&format=sam',
    headers: {'X-Remote-User': user}
  },
  function(response) {
    // Continuously update stream with data
    var body = '';
    response.on('data', function(d) {
      body += d;
    });
    response.on('end', function() {
      // Data reception is done, do whatever with it!
      console.log('RESPONSE CODE ' + response.statusCode + ' ' + response.statusMessage);
      console.log('RAW HEADERS ' + response.rawHeaders);
      console.log('RAW TRAILERS ' + response.rawTrailers);
      console.log('TRAILER data-truncated' +
        (JSON.stringify(response.trailers)));
      console.log("RECEIVED: " + body);
    });
  }
);
