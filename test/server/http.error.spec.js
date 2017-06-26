/* globals describe, it, expect, beforeAll, afterAll*/

"use strict";

const http    = require('http');
const fs      = require('fs');
const tmp     = require('tmp');

const trailer         = require('../../lib/server/http/trailer.js');
const ServerHttpError = require('../../lib/server/http/error.js');

const HttpError = ServerHttpError.HttpError;

describe('Constructor input validation', function() {
  it('response object is not given - error', function() {
    expect( () => {new HttpError();} ).toThrowError(ReferenceError,
    'HTTP response object is required');
  });
  it('response code is not given - error', function() {
    expect( () => {new HttpError({});} ).toThrowError(
    ReferenceError, 'HTTP error code is required');
  });
});

describe('Setting instance attribute values', function() {
  it('error code 500 for unknown codes', function() {
    let e = new HttpError({}, 999);
    expect(e.code).toBe(500);
    expect(e.errorType).toBe('InternalError');
    expect(e.message).toBe('Internal server error');
  });
  it('error code and type and default message', function() {
    let e = new HttpError({}, 404);
    expect(e.code).toBe(404);
    expect(e.errorType).toBe('NotFound');
    expect(e.message).toBe('Unknown error');
  });
  it('saved error message', function() {
    let e = new HttpError({}, 403, 'some error description');
    expect(e.code).toBe(403);
    expect(e.errorType).toBe('PermissionDenied');
    expect(e.message).toBe('some error description');
  });
  it('tailored status message', function() {
    let e = new HttpError({}, 400, 'some error description', 'MyCodePhrase');
    expect(e.code).toBe(400);
    expect(e.errorType).toBe('MyCodePhrase');
    expect(e.message).toBe('some error description');
  });
  it('default 400 with no rasonPhrase', function() {
    let e = new HttpError({}, 400, 'some error description');
    expect(e.code).toBe(400);
    expect(e.errorType).toBe('BadRequest');
    expect(e.message).toBe('some error description');
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

  it('Headers have been sent - set a trailer', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      trailer.declare(response);
      response.write('truncated payload');
      let e = new HttpError(response, 404);
      e. setErrorResponse();
      response.end();
    });

    http.get({socketPath: socket}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(body).toEqual('truncated payload');
        expect(response.rawTrailers).toEqual([ 'data-truncated', 'true' ]);
        done();
      });
    });
  });

  it('Headers have not been sent - set an error response', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      trailer.declare(response);
      let e = new HttpError(response, 404, 'file XX not found');
      e.setErrorResponse();
      response.end();
    });

    http.get({socketPath: socket}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(JSON.parse(body)).toEqual({
          error:   "NotFound",
          message: "file XX not found"
        });
        expect(response.statusCode).toBe(404);
        expect(response.statusMessage).toBe('file XX not found');
        let headers = response.headers;
        expect(headers['content-type']).toBe('application/json');
        expect(headers['content-length']).toBe(body.length + '');
        expect(headers['transfer-encoding']).toBe(undefined);
        expect(headers.trailer).toBe(undefined);
        expect(response.rawTrailers).toEqual([]);
        done();
      });
    });

  });

});
