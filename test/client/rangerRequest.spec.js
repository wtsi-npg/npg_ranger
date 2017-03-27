/* globals describe, expect, it, beforeAll, afterAll, jasmine */

"use strict";

const http = require('http');

const RangerRequest = require('../../lib/client/rangerRequest').RangerRequest;

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

const procjson = require('../../lib/client/rangerRequest').procJSON;

describe('JSON processing', () => {
  describe('Positives', () => {
    let json = JSON.stringify({
      format: 'BAM',
      urls:   [],
      md5:    'randommd5'
    });

    it('returns empty lists when there are no urls', () => {
      let res = procjson(json);
      expect(res.uris).toBeDefined();
      expect(res.uris.length).toEqual(0);
      expect(res.headers4uris).toBeDefined();
    });

    it('returns empty headers when no headers are passed', () => {
      let res = procjson(json);
      expect(res.headers4uris.length).toEqual(0);
    });

    it('returns individual uris for each element in json', () =>{
      let json2 = JSON.parse(json);
      json2.urls = [ {url:'url1'}, {url:'url2'}, {url:'url3'} ];
      let res = procjson(JSON.stringify(json2));
      expect(res.uris).toBeDefined();
      expect(res.uris.length).toEqual(3);
      for ( let i = 0; i < 3; i++) {
        expect(res.uris[i]).toEqual(`url${i + 1}`);
      }
    });

    it('returns an uri with headers', () => {
      let json2 = JSON.parse(json);
      json2.urls = [{
        url:'url1',
        headers: {
          "Range":         "bytes=0-1023",
          "Authorization": "Bearer xxxx"
        }
      }];
      let res = procjson(JSON.stringify(json2));
      expect(res.uris).toBeDefined();
      expect(res.uris.length).toEqual(1);
      let uri = res.uris[0];
      expect(uri).toBe('url1');
      let headers = res.headers4uris[0];
      expect(headers).toBeDefined();
      expect(headers.Range).toBeDefined();
      expect(headers.Range).toEqual('bytes=0-1023');
      expect(headers.Authorization).toBeDefined();
      expect(headers.Authorization).toEqual('Bearer xxxx');
    });
  });

  describe('Fails', () => {
    describe('Malformed urls', () => {
      it('complains about malformed urls field', () => {
        let json = JSON.stringify({
          urls: [ 'something' ]
        });
        expect(() => {
          procjson(json);
        }).toThrowError(/^Malformed JSON redirect, missing url field/);
      });
      it('complains about empty objects in urls', () => {
        let json = JSON.stringify({
          urls: [{}]
        });
        expect(() => {
          procjson(json);
        }).toThrowError(/^Malformed JSON redirect, missing url field/);
      });
      it('complains when only headers', () => {
        let json = JSON.stringify({
          urls: [{
            headers: {
              "Range":         "bytes=0-1023",
              "Authorization": "Bearer xxxx"
            }
          }]
        });
        expect(() => {
          procjson(json);
        }).toThrowError(/^Malformed JSON redirect, missing url field/);
      });
    });
  });
});
