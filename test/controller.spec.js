"use strict";

const assert  = require('assert');
const http    = require('http');
const fs      = require('fs');
const os      = require('os');
const tmp     = require('tmp');
const RangerController = require('../lib/controller.js');

describe('Creating object instance', function() {
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

  it('request object is not given - error', function() {
    expect( () => {new RangerController()} ).toThrowError(assert.AssertionError,
    'HTTP request object is required');
  });
  it('request is not an object - error', function() {
    expect( () => {new RangerController(1)} ).toThrowError(assert.AssertionError,
    'HTTP request object is required');
  });
  it('request is not an http.ClientRequest type object - error', function() {
    expect( () => {new RangerController([])} ).toThrowError(
    assert.AssertionError, 'HTTP request object is required');
  });

  let r = new http.ClientRequest();
  it('response object is not given - error', function() {
    expect( () => {new RangerController(r)} ).toThrowError(assert.AssertionError,
    'Server response object is required');
  });
  it('response is not an object - error', function() {
    expect( () => {new RangerController(r, 1)} ).toThrowError(assert.AssertionError,
    'Server response object is required');
  });
  it('response is not an http.ServerResponse type object - error', function() {
    expect( () => {new RangerController(r, [])} ).toThrowError(
    assert.AssertionError, 'Server response object is required');
  });

  it('db object is not given or is not an object - error', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect( () => {new RangerController(request, response)} ).toThrowError(
        assert.AssertionError, 'DB handle is required');
      expect( () => {new RangerController(request, response, 1)} ).toThrowError(
        assert.AssertionError, 'DB handle is required');
      response.end();
      done();      
    });
    http.get({socketPath: socket}, function(response) {});
  });

  it('Setting values of instance variables', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c;
      expect( () => {c = new RangerController(request, response, {one: "two"})} ).not.toThrow();
      expect((typeof c == 'object')).toBe(true);
      expect((c instanceof RangerController)).toBe(true);
      expect((c.request == request)).toBe(true);
      expect((c.response == response)).toBe(true);
      expect(c.db).toEqual({one: "two"});
      expect(c.tmpDir).toBe(os.tmpdir());
      expect(c.skipAuth).toBe(false);
      expect( () => {c = new RangerController(request, response, {}, 'tmp')} ).not.toThrow();
      expect(c.tmpDir).toBe('temp');
      expect(c.skipAuth).toBe(false);
      expect( () => {c = new RangerController(request, response, {}, 'tmp', 0)} ).not.toThrow();
      expect(c.tmpDir).toBe('temp');
      expect(c.skipAuth).toBe(false);
      expect( () => {c = new RangerController(request, response, {}, null, true)} ).not.toThrow();
      expect(c.tmpDir).toBe(os.tmpdir());
      expect(c.skipAuth).toBe(false);
      response.end();
      done();
    });
    http.get({socketPath: socket}, function(response) {});
  });

  it('Temporary directory should exist', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect( () => {new RangerController(request, response, {}, '/some/dir')} )
        .toThrowError(assert.AssertionError,
        "Temp data directory '/some/dir' does not exist");
    });
    http.get({socketPath: socket}, function(response) {});
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
      expect( () => {c.handleRequest()} ).toThrowError(assert.AssertionError,
        'The data host name is required');
      c.handleRequest('localhost');
      response.end();
      done();
    });
    http.get({socketPath: socket}, function(response) {});
  });

  it('Authentication error', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
    });
    http.get({socketPath: socket}, function(response) {
      var body;
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.statusCode).toEqual(401);
        expect(response.statusMessage).toEqual('Proxy authentication required');
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "InvalidAuthentication",
                   message: "Proxy authentication required"}});
        );
        done();
      });
    });
  });
});

