/* globals describe, expect, it, afterEach */

"use strict";

const http = require('http');

const RangerRequest = require('../../lib/client/rangerRequest');

describe('Testing RangerRequest public functions', () => {
  it('open parameters', ( done ) => {
    var req = new RangerRequest();

    expect(() => {req.open();}).toThrowError(/method is required/);

    var url = 'http://localhost:80/someData';

    expect(() => {req.open('GET');}).toThrowError('url is required');
    expect(() => {req.open('GET', url, 1);}).toThrowError(/async can only be boolean/);
    expect(() => {req.open('GET', url, false);}).toThrowError(/Only async requests are supported/);

    done();
  });

  it('state validation before send', ( done ) => {
    var req = new RangerRequest();
    expect(() => {req.send();}).toThrowError(/The object state must be OPENED/);
    done();
  });
});

describe('Testing RangerRequest requests', () => {
  var srv;

  afterEach(() => {
    if ( srv ) {
      srv.close();
    }
  });

  it('Client can generate a request', ( done ) => {
    var p1 = new Promise( ( resolve ) => {
      srv = http.createServer(function(request, response) {
        response.writeHead(200, {"Content-Type": "text/html"});
        response.write('<!DOCTYPE "html">');
        response.write("<html>");
        response.write("<head>");
        response.write("<title>Test result</title>");
        response.write("</head>");
        response.write("<body>");
        response.write("Test text is in the page!");
        response.write("</body>");
        response.write("</html>");
        response.end();
      });

      srv.listen(0, function() {
        resolve(srv);
      });
    });

    p1.then( ( srv ) => {
      expect(srv.address.port).not.toBe(0);

      let url = 'http://localhost:' + srv.address().port + '/someData';
      let req = new RangerRequest();
      req.open('GET', url);

      req.onreadystatechange = () => {
        if ( req.readyState === 4 ) {
          expect(req.readyState).toBe(4, 'Reached a readystate DONE');
          expect(req.statusMessage).toBe('OK', 'Correct status message');
          expect(req.status).toBe(200, 'Got a 200 request status');
          done();
        }
      };

      req.send('');
    });
  }, 3000);
});
