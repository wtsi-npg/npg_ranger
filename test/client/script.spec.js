/* globals describe, expect, it, beforeAll, afterAll, jasmine */

"use strict";

const http     = require('http');
const execFile = require('child_process').execFile;
const tmp      = require('tmp');
const fse      = require('fs-extra');

describe('Testing ranger client script', () => {

  var srv;
  var url = '';
  var page = '';
  var tmp_dir = '';

  beforeAll((done) => {

    page  = '<!DOCTYPE "html">';
    page += '<html>';
    page += '<head><title>Test result</title></head>';
    page += '<body>Test text is in the page!</body>';
    page += '</html>';

    let tmpobj = tmp.dirSync({ prefix: 'npg_ranger_test_' });
    tmp_dir = tmpobj.name;
    
    srv = http.createServer();
    srv.listen(0, function() {
      url = 'http://127.0.0.1:' + srv.address().port + '/someData';
      done();
    });
  });

  afterAll(() => {
    srv.close();
    fse.remove(tmp_dir, (err) => {
       if (err) {console.log(`Error removing ${tmp_dir}: ${err}`);}
    });
  });

  it('Client can generate a request and receive a response', ( done ) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      expect(request.headers).not.toEqual(jasmine.objectContaining({te: 'trailers'}));
      response.writeHead(200, {"Content-Type": "text/html"});
      response.write(page);
      response.end();
    });

    execFile('bin/client.js', [url], (error, stdout) => {
      expect(error).toBe(null, 'client exist normally');
      expect(stdout.toString()).toBe(page, 'response data as sent by the server');
      done();
    });

  }, 3000);

  it('Client detects server error status code and message', ( done ) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      response.writeHead(500, 'Ranger server error No 6', {"Content-Type": "text/html"});
      response.write(page);
      response.end();
    });

    execFile('bin/client.js', [url], (error, stdout) => {
      expect(error).not.toBe(null, 'client exist with an error');
      expect(error.code).toBe(1, 'error code is 1');
      expect(stdout.toString()).toBe('', 'output stream is empty');
      done();
    });

  }, 3000);

  it('Client can detect data error via a trailer', (done) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      expect(request.headers).toEqual(jasmine.objectContaining({te: 'trailers'}));
      response.setHeader('Transfer-Encoding', 'chunked');
      response.setHeader('Trailer', 'data-truncated,mytrailer,yourtrailer');
      response.writeHead(200, 'OK', {"Content-Type": "text/html"});
      response.addTrailers({'data-truncated': 'true',
                            'mytrailer':      'is a trailer',
                            'yourtrailer':    'is also a trailer'});
      response.end();
    });

    execFile('bin/client.js', [url, '--accept-trailers'], (error, stdout) => {
      expect(error).not.toBe(null, 'client exist with an error');
      expect(error.code).toBe(1, 'error code is 1');
      expect(stdout.toString()).toBe('', 'output stream is empty');
      done();
    });

  }, 3000);

  it('Client can confirm data is OK via a trailer (response empty)', (done) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      expect(request.headers).toEqual(jasmine.objectContaining({te: 'trailers'}));
      response.setHeader('Transfer-Encoding', 'chunked');
      response.setHeader('Trailer', 'data-truncated');
      response.writeHead(206, 'OK in the trailer', {"Content-Type": "text/html"});
      response.addTrailers({'data-truncated': 'false'});
      response.end();
    });

    execFile('bin/client.js', [url, '--accept-trailers'], (error, stdout) => {
      expect(error).toBe(null, 'client exist normally');
      expect(stdout.toString()).toBe('', 'output stream is empty');
      done();
    });

  }, 3000);

  it('Error if the client cannot write to the output file', ( done ) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      response.writeHead(200, {"Content-Type": "text/html"});
      response.write(page);
      response.end();
    });

    execFile('bin/client.js', [url, '/nonexisiting/partition/my.out'],
        (error) => {
      expect(error).not.toBe(null, 'client exist with an error');
      expect(error.code).toBe(1, 'error code is 1');
      done();
    });

  }, 3000);

  it('Response is written to a file', ( done ) => {
    srv.removeAllListeners('request');
    srv.on('request', (request, response) => {
      response.writeHead(200, {"Content-Type": "text/html"});
      response.write(page);
      response.end();
    });
    
    let path = tmp_dir + '/' + 'my.out';
    execFile('bin/client.js', [url, path], (error) => {
      expect(error).toBe(null, 'client exist normally');
      expect(fse.readFileSync(path).toString()).toBe(page, 'output file content is correct');
      done();
    });

  }, 3000);

});
