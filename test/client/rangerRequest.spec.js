/* globals describe, expect, it, beforeAll, afterAll, jasmine */

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
  var url = '';
  var page = '';

  beforeAll((done) => {

    page  = '<!DOCTYPE "html">';
    page += '<html>';
    page += '<head><title>Test result</title></head>';
    page += '<body>Test text is in the page!</body>';
    page += '</html>';

    srv = http.createServer();
    srv.listen(0, function() {
      url = 'http://127.0.0.1:' + srv.address().port + '/someData';
      done();
    });
  });

  afterAll(() => {
    srv.close();
  });

  it('Client can generate a request and receive a response', ( done ) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      // construct RangerRequest with acceptTrailers = false, so TE: trailers header
      // should not exist
      expect(request.headers).not.toEqual(jasmine.objectContaining({te: 'trailers'}));
      response.writeHead(200, {"Content-Type": "text/html"});
      response.write(page);
      response.end();
    });

    let req = new RangerRequest();
    req.open('GET', url);
    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.statusMessage).toBe('OK', 'Correct status message');
        expect(req.status).toBe(200, 'Got a 200 request status');
        expect(req.response.toString()).toBe(page, 'response data as sent by the server');
        done();
      }
    };

    req.send('');

  }, 3000);

  it('Client detects server error status code and message', ( done ) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      expect(request.headers).not.toEqual(jasmine.objectContaining({te: 'trailers'}));
      response.writeHead(500, 'Ranger server error No 6', {"Content-Type": "text/html"});
      response.write(page);
      response.end();
    });

    let req = new RangerRequest();
    req.open('GET', url);
    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.statusMessage).toBe(
            'Ranger server error No 6 for ' + url, 'Correct status message');
        expect(req.status).toBe(500, 'Request status 500 as set by the server');
        expect(req.response).toBe(undefined, 'No data in response');
        done();
      }
    };

    req.send();

  }, 3000);

  it('Client can detect data error via a trailer', (done) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      expect(request.headers).toEqual(jasmine.objectContaining({te: 'trailers'}));
      response.setHeader('Transfer-Encoding', 'chunked');
      response.setHeader('Trailer', 'data-truncated,mytrailer,yourtrailer');
      response.writeHead(200, 'OK', {"Content-Type": "text/html"});
      response.addTrailers({'data-truncated': 'true',
                            mytrailer:        'is a trailer',
                            yourtrailer:      'is also a trailer'});
      response.end();
    });

    let req = new RangerRequest(true);
    req.open('GET', url);
    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.statusMessage).toBe(
          'Incomplete or truncated data for ' + url, 'Data truncation status message');
        expect(req.status).toBe(424, 'The client changes status code from 200 to 424');
        done();
      }
    };

    req.send('');

  }, 3000);

  it('Client can confirm data is OK via a trailer', (done) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      expect(request.headers).toEqual(jasmine.objectContaining({te: 'trailers'}));
      response.setHeader('Transfer-Encoding', 'chunked');
      response.setHeader('Trailer', 'data-truncated');
      response.writeHead(206, 'OK in the trailer', {"Content-Type": "text/html"});
      response.write(page);
      response.addTrailers({'data-truncated': 'false'});
      response.end();
    });

    let req = new RangerRequest(true);
    req.open('GET', url);
    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.statusMessage).toBe('OK in the trailer', 'Status message');
        expect(req.status).toBe(206, 'Status code 206');
        expect(req.response.toString()).toBe(page, 'response as sent by the server');
        done();
      }
    };

    req.send('');

  }, 3000);

  it('If not accepting trailers, than not checking trailers either', (done) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      response.setHeader('Transfer-Encoding', 'chunked');
      response.setHeader('Trailer', 'data-truncated');
      response.writeHead(200, 'Fail in the trailer',{"Content-Type": "text/html"});
      response.write(page);
      let h = {};
      h['data-truncated'] = 'true';
      response.addTrailers(h);
      response.end();
    });

    let req = new RangerRequest();
    req.open('GET', url);
    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.statusMessage).toBe('Fail in the trailer', 'Status message');
        expect(req.status).toBe(200, 'Status code 200');
        expect(req.response.toString()).toBe(page, 'response as sent by the server');
        done();
      }
    };

    req.send();

  }, 3000);

  it('Not getting trailers when asked is not a failure', (done) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      expect(request.headers).toEqual(jasmine.objectContaining({te: 'trailers'}));
      response.setHeader('Transfer-Encoding', 'chunked');
      response.writeHead(200, 'OK, no trailer', {"Content-Type": "text/html"});
      response.write(page);
      response.end();
    });

    let req = new RangerRequest(true);
    req.open('GET', url);
    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.statusMessage).toBe('OK, no trailer', 'Status message');
        expect(req.status).toBe(200, 'Status code 200');
        expect(req.response.toString()).toBe(page, 'response as sent by the server');
        done();
      }
    };

    req.send();

  }, 3000);
});
