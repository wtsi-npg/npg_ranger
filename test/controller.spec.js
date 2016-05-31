/* globals describe, it, expect, beforeAll, afterAll*/

"use strict";

const assert  = require('assert');
const http    = require('http');
const fs      = require('fs');
const os      = require('os');
const tmp     = require('tmp');
const RangerController = require('../lib/controller.js');

describe('Creating object instance - synch', function() {
  it('request object is not given - error', function() {
    expect( () => {new RangerController();} ).toThrowError(assert.AssertionError,
    'HTTP request object is required');
  });
  it('request is not an object - error', function() {
    expect( () => {new RangerController(1);} ).toThrowError(assert.AssertionError,
    'HTTP request object is required');
  });
  it('response object is not given - error', function() {
    expect( () => {new RangerController({});} ).toThrowError(assert.AssertionError,
    'Server response object is required');
  });
  it('response is not an object - error', function() {
    expect( () => {new RangerController({}, 1);} ).toThrowError(assert.AssertionError,
    'Server response object is required');
  });
  it('response is not an http.ServerResponse type object - error', function() {
    expect( () => {new RangerController({}, []);} ).toThrowError(
    assert.AssertionError, 'Server response object is required');
  });
});

describe('set error response', function() {
  // Create server HTTP server object.
  const server = http.createServer();
  // Generate synchronously a temporary file name.
  var socket = tmp.tmpNameSync();

  beforeAll(function() {
    // Start listening on a socket
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
    });
  });

  // This tidy-up callback is not called when the spec exits early
  // due to an error. Seems to be a bug in jasmine.
  afterAll(function() {
    server.close();
    try { fs.unlinkSync(socket); } catch (e) {}
  });

  it('db object is not given or is not an object - error', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      console.log(response.toString());
      assert(typeof request == 'object');
      expect( () => {new RangerController(request, response);} ).toThrowError(
        assert.AssertionError, 'DB handle object is required');
      expect( () => {new RangerController(request, response, 1);} ).toThrowError(
        assert.AssertionError, 'DB handle object is required');
      response.write('payload');
      response.end();
    });

    http.get({socketPath: socket}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        done();
      });
    });
  });

  it('Setting values of instance variables', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c;
      expect( () => {c = new RangerController(request, response, {one: "two"});} ).not.toThrow();
      expect((typeof c == 'object')).toBe(true);
      expect((c instanceof RangerController)).toBe(true);
      expect((c.request == request)).toBe(true);
      expect((c.response == response)).toBe(true);
      expect(c.db).toEqual({one: "two"});
      expect(c.tmpDir).toBe(os.tmpdir());
      expect(c.skipAuth).toBe(false);
      expect( () => {c = new RangerController(request, response, {}, null, 0);} ).not.toThrow();
      expect(c.tmpDir).toBe(os.tmpdir());
      expect(c.skipAuth).toBe(false);
      expect( () => {c = new RangerController(request, response, {}, '', true);} ).not.toThrow();
      expect(c.tmpDir).toBe(os.tmpdir());
      expect(c.skipAuth).toBe(true);
      response.end();
      done();
    });
    http.get({socketPath: socket}, function() {});
  });

  it('Temporary directory should exist', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect( () => {new RangerController(request, response, {}, '/some/dir');} )
        .toThrowError(assert.AssertionError,
        "Temp data directory '/some/dir' does not exist");
      let c;
      expect( () => {c = new RangerController(request, response, {}, 'test', 0);} ).not.toThrow();
      expect(c.tmpDir).toBe('test');
      done();
    });
    http.get({socketPath: socket}, function() {});
  });
});

describe('Handling requests - error responses', function() {
  const server = http.createServer();
  var socket = tmp.tmpNameSync();
  beforeAll(function() {
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
    });
  });
  afterAll(function() {
    server.close();
    try { fs.unlinkSync(socket); } catch (e) {}
  });

  it('Data host name argument is required', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).toThrowError(assert.AssertionError,
        'The data host name is required');
      c.handleRequest('localhost');
      response.end();
      done();
    });
    http.get({socketPath: socket}, function() {});
  });

  it('Authentication error', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect(c.skipAuth).toBe(false);
      expect( () => {c.handleRequest('localhost');} ).not.toThrow();
    });
    http.get({socketPath: socket}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(401);
        expect(response.statusMessage).toEqual('Proxy authentication required');
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "InvalidAuthentication",
                   message: "Proxy authentication required"}});
        done();
      });
    });
  });

  it('Not found error, no auth', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"}, null, true);
      expect(c.skipAuth).toBe(true);
      expect( () => {c.handleRequest('localhost');} ).not.toThrow();
    });
    http.get({socketPath: socket, path: '/invalid'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : /invalid');
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "NotFound",
                   message: "URL not found : /invalid"}});
        done();
      });
    });
  });

  it('Not found error, auth checked', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"}, null, false);
      expect(c.skipAuth).toBe(false);
      expect( () => {c.handleRequest('localhost');} ).not.toThrow();
    });
    let req = http.request({socketPath: socket, path: '/invalid'});
    req.setHeader('X-Remote-User', 'user1');
    req.end();
    req.on('response', function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : /invalid');
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "NotFound",
                   message: "URL not found : /invalid"}});
        done();
      });
    });
  });

  it('Invalid input error for a sample url', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"}, null, true);
      expect(c.skipAuth).toBe(true);
      expect( () => {c.handleRequest('localhost');} ).not.toThrow();
    });

    http.get({socketPath: socket, path: '/sample'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        let m = 'Invalid request: sample accession number should be given';
        expect(response.statusMessage).toEqual(m);
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "InvalidInput",
                   message: m}});
        done();
      });
    });
  });

  it('Invalid input error for a file url', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"}, null, true);
      expect(c.skipAuth).toBe(true);
      expect( () => {c.handleRequest('localhost');} ).not.toThrow();
    });

    http.get({socketPath: socket, path: '/file'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        let m = 'Invalid request: file name should be given';
        expect(response.statusMessage).toEqual(m);
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "InvalidInput",
                   message: m}});
        done();
      });
    });
  });
});