/* globals describe, it, expect, beforeAll, afterAll*/

"use strict";

const http    = require('http');
const fs      = require('fs');
const tmp     = require('tmp');
const trailer = require('../../lib/server/http/trailer.js');

describe('Input validation', function() {
  it('declare - response object is not given - error', function() {
    expect( () => {trailer.declare();} ).toThrowError(ReferenceError,
    'HTTP response object is required');
  });
  it('setDataTruncation - response object is not given - error', function() {
    expect( () => {trailer.setDataTruncation();} ).toThrowError(
    ReferenceError, 'HTTP response object is required');
  });
  it('setDataTruncation - response object is not given - error', function() {
    expect( () => {trailer.setDataTruncation({});} ).toThrowError(
    ReferenceError, 'boolean flag indicating data truncation is required');
  });
  it('setDataTruncation - response object is not given - error', function() {
    expect( () => {trailer.setDataTruncation({}, 4);} ).toThrowError(
    ReferenceError, 'boolean flag indicating data truncation is required');
  });
});

describe('declaring, setting and removing a trailer', function() {
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

  it('Declare and set a trailer to mark data truncation', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      trailer.declare(response);
      expect(response.getHeader('Trailer')).toBe('data-truncated');
      response.write('useful payload');
      trailer.setDataTruncation(response, true);
      response.end();
    });

    http.get({socketPath: socket}, function(response) {
      response.on('data', function() {
        // not interested in data, but the end event is not called
        // unless the data is processed
      });
      response.on('end', function() {
        expect(response.rawTrailers).toEqual([ 'data-truncated', 'true' ]);
        done();
      });
    });
  });

  it('Declare and set a trailer to mark good data', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      trailer.declare(response);
      expect(response.getHeader('Trailer')).toBe('data-truncated');
      response.write('useful payload');
      trailer.setDataTruncation(response, false);
      response.end();
    });

    http.get({socketPath: socket}, function(response) {
      response.on('data', function() {});
      response.on('end', function() {
        expect(response.rawTrailers).toEqual([ 'data-truncated', 'false' ]);
        done();
      });
    });

  });

  it('Declare and remove a trailer', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      trailer.declare(response);
      expect(response.getHeader('Trailer')).toBe('data-truncated');
      expect( () => {trailer.removeDeclaration(response);} ).not.toThrow();
      expect(response.getHeader('Trailer')).toBe(undefined);
      response.end();
      done();
    });

    http.get({socketPath: socket}, function(response) {
      response.on('data', function() {});
      response.on('end', function() {});
    });

  });

  it('Trailer has not been declared: no error removing declaration, error setting', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect( () => {trailer.removeDeclaration(response);} ).not.toThrow();
      response.write('useful payload');
      expect( () => {trailer.setDataTruncation(response, true);} )
        .toThrowError(Error,
        'Cannot set data truncation trailer because it has not been declared');
      response.end();
    });

    http.get({socketPath: socket}, function(response) {
      response.on('data', function() {});
      response.on('end', function() {
        expect(response.rawTrailers).toEqual([]);
        done();
      });
    });

  });

  it('If Transfer-Encoding header is not set, Transfer-Encoding is set to "chunked"', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      response.removeHeader('Transfer-Encoding');
      response.setHeader('Content-Type', 'application/json');
      trailer.declare(response);
      expect(response.getHeader('Transfer-Encoding')).toBe('chunked');
      expect(response.getHeader('Trailer')).toBe('data-truncated');
      response.write('{"some": "property"}');
      trailer.setDataTruncation(response, true);
      response.end();
    });

    http.get({socketPath: socket}, function(response) {
      response.on('data', function() {});
      response.on('end', function() {
        let rawTrailers = response.rawTrailers;
        expect(rawTrailers.length).toEqual(2);
        expect(rawTrailers[0]).toEqual('data-truncated');
        expect(rawTrailers[1]).toEqual('true');
        done();
      });
    });

  });

  it('Trailer cannot be declared after the headers have been sent', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      response.write('payload');
      expect( () => {trailer.declare(response);} ).toThrowError(Error,
        "Cannot set headers after they are sent to the client");
      expect(response.getHeader('Trailer')).toBe(undefined);
      response.end();
      done();
    });

    http.get({socketPath: socket}, function(response) {
      response.on('data', function() {});
      response.on('end',  function() {});
    });

  });

});
