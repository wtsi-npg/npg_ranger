/* globals describe, it, xit, expect, beforeAll, beforeEach, afterAll, jasmine */

"use strict";

const assert  = require('assert');
const http    = require('http');
const fs      = require('fs');
const fse     = require('fs-extra');
const tmp     = require('tmp');
const RangerController = require('../../lib/server/controller.js');
const config  = require('../../lib/config.js');
const constants = require('../../lib/constants.js');
const tokenUtils = require('../../lib/token_utils.js');

const utils           = require('./test_utils.js');
const ServerHttpError = require('../../lib/server/http/error');
const trailer         = require('../../lib/server/http/trailer.js');

const GA4GH_URL       = '/ga4gh/sample';
const GA4GH_TOKEN_URL = '/authtoken' + GA4GH_URL;

// Create temp dir here so it is available for all tests.
// Use this dir as a default dir that will be available in all.
var tmpDir    = config.tempFilePath('npg_ranger_controller_test_');
var dummy     = function() { return {tempdir: tmpDir, skipauth: true}; };
var options;

let checkGenericErrorResponse = ( expect, response ) => {
  let headers = response.headers;
  expect(headers['content-type']).toEqual('application/json');
  expect(headers['transfer-encoding']).toBe(undefined);
  expect(headers.trailer).toBe(undefined);
  expect(response.rawTrailers).toEqual([]);
};

describe('Creating object instance - synch', function() {
  beforeAll(function() {
    options = config.provide(dummy);
    fse.ensureDirSync(tmpDir);
  });

  afterAll( () => {
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });

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
});

describe('set data truncation', () => {
  const server = http.createServer();
  var socket = tmp.tmpNameSync();

  beforeAll( ( done ) => {
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });

  afterAll( () => {
    server.close();
    utils.removeSocket(socket);
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });

  beforeEach( () => {
    server.removeAllListeners('request');
  });

  let requestOptions = [
    {method: 'GET', socketPath: socket},
    {method: 'POST', socketPath: socket}
  ];
  
  requestOptions.forEach(obj => {
    it(`sets sends error response when headers have not been sent for ${obj.method}`, ( done ) => {
      server.on('request', ( request, response ) => {
        assert(typeof request === 'object');
        trailer.declare(response);
        let c = new RangerController(request, response, {one: "two"});
        c.errorResponse(response, 404, 'Not Found');
      });

      http.request( obj, ( response ) => {
        var body = '';
        response.on('data', d => { body += d;});
        response.on('end', () => {
          expect(response.statusCode).toEqual(404);
          let jsonbody = JSON.parse(body);
          expect(jsonbody.htsget).toBeDefined();
          let htsget = jsonbody.htsget;
          expect(htsget.error).toEqual('NotFound');
          expect(htsget.message).toEqual('Not Found');
          expect(response.rawTrailers).toEqual([]);
          done();
        });
      }).end();
    });

    it(`sets up trailers when setting error response after headers sent for ${obj.method}`, ( done ) => {
      server.on('request', ( request, response ) => {
        assert(typeof request == 'object');
        trailer.declare(response);
        let c = new RangerController(request, response, {one: "two"});
        response.write('truncated payload'); //send headers
        c.errorResponse(response, 404);
      });

      http.request( obj, ( response ) => {
        var body = '';
        response.on('data', d => { body += d;});
        response.on('end', () => {
          expect(response.statusCode).toEqual(200);
          expect(body).toEqual('truncated payload');
          expect(response.rawTrailers).toEqual([ 'data-truncated', 'true' ]);
          done();
        });
      }).end();
    });
  });
});

describe('set error response', function() {
  // Create server HTTP server object.
  const server = http.createServer();
  // Generate synchronously a temporary file name.
  var socket = tmp.tmpNameSync();

  beforeAll((done) => {
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
    // Start listening on a socket
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });

  // This tidy-up callback is not called when the spec exits early
  // due to an error. Seems to be a bug in jasmine.
  afterAll( () => {
    server.close();
    utils.removeSocket(socket);
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });

  let requestOptions = [
    {method: 'GET', socketPath: socket},
    {method: 'POST', socketPath: socket}
  ];

  requestOptions.forEach( obj => {
    it(`db object is not given or is not an object - error for ${obj.method}`, function(done) {
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        assert(typeof request == 'object');
        expect( () => {new RangerController(request, response);} ).toThrowError(
          assert.AssertionError, 'DB handle object is required');
        expect( () => {new RangerController(request, response, 1);} ).toThrowError(
          assert.AssertionError, 'DB handle object is required');
        response.write('payload');
        response.end();
      });
      
      http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          done();
        });
      }).end();
    });

    it(`Setting values of instance variables for ${obj.method}`, function(done) {
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        let c;
        expect( () => {c = new RangerController(request, response, {one: "two"});} ).not.toThrow();
        expect((typeof c == 'object')).toBe(true);
        expect((c instanceof RangerController)).toBe(true);
        expect((c.request == request)).toBe(true);
        expect((c.response == response)).toBe(true);
        expect(c.db).toEqual({one: "two"});
        response.end();
        done();
      });
      http.request( obj, function() {}).end();
    });
  });
});

describe('Handling requests - error responses', function() {
  const server = http.createServer();
  var socket = tmp.tmpNameSync();
  beforeAll((done) => {
    options = config.provide(dummy);
    fse.ensureDirSync(tmpDir);
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });

  afterAll( () => {
    server.close();
    try {
      fs.access(socket, fs.constants.W_OK, ( err ) => {
        if ( !err ) {
          fs.unlinkSync(socket);
        } else {
          console.log( err );
        }
      });
    } catch (e) { console.log(e); }
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });

  let requestOptions = [
    {method: 'GET', socketPath: socket, path: ''},
    {method: 'POST', socketPath: socket, path: ''}
  ];


  xit('Method not allowed error', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });

    let options = { socketPath: socket, path:'/file', method: 'POST' };
    let req = http.request(options);
    req.on('response', (response) => {
      checkGenericErrorResponse(expect, response);
      expect(response.statusCode).toEqual(400);
      done();
    });
    req.end();
  });

  requestOptions.forEach(obj => {
    it(`Authentication error for ${obj.method}`, function(done) {
      config.provide( () => { 
        return {tempdir:  tmpDir,
                skipauth: false};
      });
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        let c = new RangerController(request, response, {one: "two"});
        expect( () => {c.handleRequest();} ).not.toThrow();
      });
      obj.path = '/authuser';
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(401);
          expect(response.statusMessage).toEqual('Proxy authentication required');
          expect(JSON.parse(body)).toEqual({
            htsget: {
              error:   "InvalidAuthentication",
              message: "Proxy authentication required"
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`Not found error, no auth for ${obj.method}`, function(done) {
      config.provide( () => {
        return {tempdir:  tmpDir,
                skipauth: true};
      });
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        let c = new RangerController(request, response, {one: "two"});
        expect( () => {c.handleRequest();} ).not.toThrow();
      });
      obj.path = '/invalid';
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(404);
          expect(response.statusMessage).toEqual('URL not found : /invalid');
          expect(JSON.parse(body)).toEqual({
            htsget:{
              error:   "NotFound",
              message: "URL not found : /invalid"
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`Not found error, auth checked for ${obj.method}`, function(done) {
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        let c = new RangerController(request, response, {one: "two"});
        expect( () => {c.handleRequest();} ).not.toThrow();
      });
      obj.path = '/invalid';
      let req = http.request( obj );
      req.setHeader('X-Remote-User', 'user1');
      req.on('response', function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(404);
          expect(response.statusMessage).toEqual('URL not found : /invalid');
          expect(JSON.parse(body)).toEqual({
            htsget: {
              error:   "NotFound",
              message: "URL not found : /invalid"
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`Invalid input error for duplicate values for an attribute for ${obj.method}`, ( done ) => {
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        let c = new RangerController(request, response, {one: "two"});
        expect( () => {c.handleRequest();} ).not.toThrow();
      });
      if (obj.method === "POST") {
        obj.path = '/sample';
      } else {
        obj.path = '/sample?attr1=value1&attr2=value2&attr1=value3';
      }
      let req = http.request( obj, ( response ) => {
        var body = '';
        response.on('data', ( d ) => { body += d;});
        response.on('end', () => {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(400);
          let m = "Invalid request: multiple values for attribute 'attr1'";
          expect(response.statusMessage).toEqual(m);
          expect(JSON.parse(body)).toEqual({
            htsget: {
              error:   ServerHttpError.INVALID_INPUT,
              message: m
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write('{"format":"bam","attr1":"value","attr2":"value2","attr1":"value3"}', () => {
          req.end();
        });
      } else {
        req.end();
      }
    });


    it(`Invalid input error for a sample url for ${obj.method}`, ( done ) => {
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        let c = new RangerController(request, response, {one: "two"});
        expect( () => {c.handleRequest();} ).not.toThrow();
      });
      obj.path = '/sample';
      let req = http.request( obj, ( response ) => {
        var body = '';
        response.on('data', ( d ) => { body += d;});
        response.on('end', () => {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(400);
          let m = 'Invalid request: sample accession number should be given';
          expect(response.statusMessage).toEqual(m);
          expect(JSON.parse(body)).toEqual({
            htsget: {
              error:   ServerHttpError.INVALID_INPUT,
              message: m
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`Invalid input error for a file url for ${obj.method}`, ( done ) => {
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        let c = new RangerController(request, response, {one: "two"});
        expect( () => {c.handleRequest();} ).not.toThrow();
      });
      obj.path = '/file';
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(400);
          let m = 'Invalid request: file name should be given';
          expect(response.statusMessage).toEqual(m);
          expect(JSON.parse(body)).toEqual({
            htsget: {
              error:   ServerHttpError.INVALID_INPUT,
              message: m
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`Invalid input error for a vcf file when multiref set for ${obj.method}`, (done) => {

      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        // Set multiref mode
        config.provide(() => { return {tempdir: tmpDir, multiref: true, skipauth: true}; });
        let c = new RangerController(request, response, {one: "two"});
        expect( () => {c.handleRequest();} ).not.toThrow();
        // unset multiref
        config.provide(dummy);
      });
      if (obj.method === "POST") {
        obj.path = '/sample';
      } else {
        obj.path = '/sample?accession=XYZ120923&format=vcf';
      }
      let req = http.request( obj, (response) => {
        var body = '';
        response.on('data', (d) => { body += d;});
        response.on('end', () => {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(400);
          let m = 'Invalid request: cannot produce VCF files while multiref set on server';
          expect(response.statusMessage).toEqual(m);
          expect(JSON.parse(body)).toEqual({
            htsget: {
              error:   ServerHttpError.INVALID_INPUT,
              message: m
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"vcf"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });
  });
});

describe('Handling POST requests', function() {
  const server = http.createServer();
  var socket = tmp.tmpNameSync();
  beforeAll((done) => {
    options = config.provide(dummy);
    fse.ensureDirSync(tmpDir);
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });

  afterAll( () => {
    server.close();
    try {
      fs.access(socket, fs.constants.W_OK, ( err ) => {
        if ( !err ) {
          fs.unlinkSync(socket);
        } else {
          console.log( err );
        }
      });
    } catch (e) { console.log(e); }
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });

  it('POST connection valid', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      request.on('data', (chunk) => {
        console.log(JSON.parse(chunk));
      });
      let c = new RangerController(request, response, {});
      expect( () => {c.handleRequest();} ).not.toThrow();
      // Creates the server
    });

    let options = {
      socketPath: socket,
      path:"/file",
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      }
    };

    let testVar = {
      "format":"bam"
    };

    let req = http.request(options);
    // create client to connect afterwards
    req.write(JSON.stringify(testVar));
    req.on('response', (response) => {
      checkGenericErrorResponse(expect, response);
      expect(response.statusCode).toEqual(400);
      expect(response.statusMessage).toMatch('file name should be given');
      done();
    });
    req.end();
  });
});

describe('Sample reference', () => {
  const server = http.createServer();
  let   socket = tmp.tmpNameSync();
  let   id     = 'XYZ120923';

  const child       = require('child_process');
  const MongoClient = require('mongodb').MongoClient;

  const BASE_PORT  = 1400;
  const PORT_RANGE = 200;
  const PORT       = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;
  const FIXTURES   = 'test/server/data/fixtures/fileinfo.json';

  afterAll( () => {
    server.close();

    child.execSync(`mongo 'mongodb://localhost:${PORT}/admin' --eval 'db.shutdownServer()'`);
    console.log('\nMONGODB server has been shut down');

    utils.removeSocket(socket);
  });

  beforeAll( (done) => {
    options = config.provide(dummy);
    fse.ensureDirSync(tmpDir);

    console.log(`MONGO data directory: ${tmpDir}`);
    let db_name = 'npg_ranger_test';
    let url = `mongodb://localhost:${PORT}/${db_name}`;

    let command = `mongod -f test/server/data/mongodb_conf.yml --port ${PORT} --dbpath ${tmpDir} --pidfilepath ${tmpDir}/mpid --logpath ${tmpDir}/dbserver.log`;
    console.log(`\nCommand to start MONGO DB daemon: ${command}`);
    let out = child.execSync(command);
    console.log(`Started MONGO DB daemon: ${out}`);
    command = `mongoimport --port ${PORT} --db ${db_name} --collection fileinfo --jsonArray --file ${FIXTURES}`;
    out = child.execSync(command);
    console.log(`Loaded data to MONGO DB: ${out}`);

    MongoClient.connect(url, (err, db) => {
      assert.equal(err, null);
      server.on('request', (request, response) => {
        let c = new RangerController(request, response, db);
        c.handleRequest();
      });
      server.listen(socket, () => {
        console.log(`Server listening on socket ${socket}`);
        done();
      });
    });
  });

  [
    '/sample//reference',
    '/sample/-/reference',
    '/sample/%20/reference', // Encoded blank space
    `/sample/${id}/reference/something`
  ].forEach( ( thisPath ) => {
    it(`returns error response for invalid url '${thisPath}'`, ( done ) => {
      http.get({
        socketPath: socket,
        path:       thisPath
      }, ( response ) => {
        let body = '';
        response.on('data', ( d ) => { body += d;});
        response.on('end',  ()    => {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(404);
          expect(response.statusMessage).toEqual(`URL not found : ${thisPath}`);
          done();
        });
      });
    });
  });

  it('gets reference for existing accession', ( done ) => {
    let thisPath = `/sample/${id}/reference`;
    http.get({
      socketPath: socket,
      path:       thisPath
    }, ( response ) => {
      let body = '';
      response.on('data', ( d ) => {
        body += d;
      });
      response.on('end',  ()    => {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(`OK`);
        expect(JSON.parse(body)).toEqual(
          jasmine.objectContaining({
            reference: '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa',
            accession: 'XYZ120923'
          })
        );
        done();
      });
    });
  });

  it('gets 404 when accession does not exist', ( done ) => {
    let missingAcc = 'NONE';
    let thisPath   = `/sample/${missingAcc}/reference`;
    http.get({
      socketPath: socket,
      path:       thisPath
    }, ( response ) => {
      let body = '';
      response.on('data', ( d ) => { body += d;});
      response.on('end',  ()    => {
        checkGenericErrorResponse(expect, response);
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual(
          `No files for sample accession ${missingAcc}`
        );
        done();
      });
    });
  });
});

describe('Redirection with token in json response', () => {
  const server = http.createServer();
  var socket   = tmp.tmpNameSync();
  let id       = 'EGA45678';
  let server_path = `${GA4GH_TOKEN_URL}/${id}`;

  beforeAll((done) =>  {
    fse.ensureDirSync(tmpDir);
    options = config.provide(() => {
      return { tempdir: tmpDir };
    });
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {});
      c.handleRequest();
    });
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });

  afterAll( () => {
    server.close();
    utils.removeSocket(socket);
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });

  let requestOptions = [
    {method: 'GET', socketPath: socket, path: server_path},
    {method: 'POST', socketPath: socket, path: server_path}
  ];
  
  requestOptions.forEach( obj => {
    it('successful redirection with token', function(done) {
      let my_token = 'XXXYYYXXX';
      let headers  = {};
      headers[
        constants.TOKEN_BEARER_KEY_NAME
      ] = tokenUtils.formatTokenForHeader(my_token);
      obj.headers = headers;
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          expect(response.headers['content-type']).toMatch(
            /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
          );
          expect(response.statusCode).toEqual(200);
          expect(response.statusMessage).toEqual(
            'OK, see redirection instructions in the body of the message');
          let url = `http://localhost/sample?accession=${id}&format=BAM`;
          expect(JSON.parse(body)).toEqual({
            htsget: {
              format: 'BAM',
              urls: [
                {
                  'url': url,
                  'headers': {
                    'Authorization': `Bearer ${my_token}`
                  }
                }
              ]
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });
  });
});

describe('Redirection in json response', function() {
  const server = http.createServer();
  var socket   = tmp.tmpNameSync();
  let id       = 'EGA45678';
  let server_path = `${GA4GH_URL}/${id}`;

  beforeAll((done) =>  {
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {});
      c.handleRequest();
    });
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });

  afterAll( () => {
    server.close();
    utils.removeSocket(socket);
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });


  it(`successful redirection, filter given`, function(done) {
    options = {method: "GET", socketPath: socket, path : server_path + '?target=0&manual_qc=&alignment_not=undef'};
    let req = http.request( options, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toMatch(
          /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
        );
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM&target=0&manual_qc=&alignment_not=undef`;
        expect(JSON.parse(body)).toEqual({
          htsget: {
            format: 'BAM',
            urls: [{'url': url}]
          }
        });
        done();
      });
    });
    req.end();
  });
        
  let requestOptions = [
    {method: 'GET', socketPath: socket},
    {method: 'POST', socketPath: socket}
  ];

  it('redirection error, invalid character in reference name for GET', function(done) {
    let obj = {method: 'GET', socketPath: socket, path: server_path + '?referenceName=chr@1&start=4&end=400'};
    let req = http.request ( obj, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        checkGenericErrorResponse(expect, response);
        expect(response.statusCode).toEqual(400);
        let errorPayload = JSON.parse(body);
        expect(errorPayload.htsget).toBeDefined();
        expect(errorPayload.htsget.error).toBe(ServerHttpError.INVALID_INPUT);
        expect(response.statusMessage).toEqual(
          'Invalid character in reference name chr@1');
        done();
      });
    });
    req.end();
  });
  
  requestOptions.forEach( optionsList => {
    it(`invalid url - no id - error response for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      obj.path = GA4GH_URL;
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(404);
          expect(response.statusMessage).toEqual('URL not found : ' + GA4GH_URL);
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });
    
    it(`invalid url - no id - error response for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      let path = GA4GH_URL + '/';
      obj.path = path;
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(404);
          expect(response.statusMessage).toEqual('URL not found : ' + path);
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`invalid sample id - error response for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      let path = GA4GH_URL + 'ERS-4556';
      obj.path = path;
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(404);
          expect(response.statusMessage).toEqual('URL not found : ' + path);
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });
    
    it(`successful redirection, no query params for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      obj.path = server_path;
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          expect(response.headers['content-type']).toMatch(
            /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
          );
          expect(response.statusCode).toEqual(200);
          expect(response.statusMessage).toEqual(
            'OK, see redirection instructions in the body of the message');
          let url = `http://localhost/sample?accession=${id}&format=BAM`;
          expect(JSON.parse(body)).toEqual({
            htsget: {
              format: 'BAM',
              urls: [{'url': url}]
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`successful redirection, format given for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      obj.path = server_path + '?format=CRAM';
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          expect(response.headers['content-type']).toMatch(
            /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
          );
          expect(response.statusCode).toEqual(200);
          expect(response.statusMessage).toEqual(
            'OK, see redirection instructions in the body of the message');
          let url = `http://localhost/sample?accession=${id}&format=CRAM`;
          expect(JSON.parse(body)).toEqual({
            htsget: {
              format: 'CRAM',
              urls: [{'url': url}]
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"cram"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`successful redirection, chromosome given for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      if (obj.method === "GET") {
        obj.path = server_path + '?referenceName=chr1';
      } else {
        obj.path = server_path;
      }
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          expect(response.headers['content-type']).toMatch(
            /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
          );
          expect(response.statusCode).toEqual(200);
          expect(response.statusMessage).toEqual(
            'OK, see redirection instructions in the body of the message');
          if (obj.method === "GET") {
            let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1`;
            expect(JSON.parse(body)).toEqual({
              htsget: {
                format: 'BAM',
                urls: [{'url': url}]
              }
            });
          } else {
            let url = `http://localhost/sample?accession=${id}&format=BAM`;
            let values = JSON.parse(body);
            expect(values.htsget).toBeDefined();
            expect(values.htsget.urls).toBeDefined();
            expect(values.htsget.format).toEqual("BAM");
            expect(values.htsget.urls[0].url).toEqual(url);
            expect(values.htsget.urls[0].headers).toEqual({"encoded_regions": "H4sIAAAAAAAAA4uuVipKTUstSs1LTvVLzE1VslJKzigyVKqNBQDCypWAGgAAAA=="});
          }
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam", "regions":[{"referenceName":"chr1"}]}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`successful redirection, range start given for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      if (obj.method === "GET") {
        obj.path = server_path + '?referenceName=chr1&start=3';
      } else {
        obj.path = server_path;
      }
      let req = http.request( obj , function(response) { 
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          expect(response.headers['content-type']).toMatch(
            /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
          );
          expect(response.statusCode).toEqual(200);
          expect(response.statusMessage).toEqual(
            'OK, see redirection instructions in the body of the message');
          if (obj.method === "GET") {
            let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1%3A4`;
            expect(JSON.parse(body)).toEqual({
              htsget: {
                format: 'BAM',
                urls: [{'url': url}]
              }
            });
          } else {
            let url = `http://localhost/sample?accession=${id}&format=BAM`;
            let values = JSON.parse(body);
            expect(values.htsget).toBeDefined();
            expect(values.htsget.urls).toBeDefined();
            expect(values.htsget.format).toEqual("BAM");
            expect(values.htsget.urls[0].url).toEqual(url);
            expect(values.htsget.urls[0].headers).toEqual({"encoded_regions": "H4sIAAAAAAAAA4uuVipKTUstSs1LTvVLzE1VslJKzigyVNJRKi5JLCpRsjKujQUAJUieXSQAAAA="});
          }
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam", "regions":[{"referenceName":"chr1", "start":3}]}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`successful redirection, range end given for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      if (obj.method === "GET") {
        obj.path = server_path + '?referenceName=chr1&end=4';
      } else {
        obj.path = server_path;
      }
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          expect(response.headers['content-type']).toMatch(
            /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
          );
          expect(response.statusCode).toEqual(200);
          expect(response.statusMessage).toEqual(
            'OK, see redirection instructions in the body of the message');
          if (obj.method === "GET") {
            let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1%3A1-4`;
            expect(JSON.parse(body)).toEqual({
              htsget: {
                format: 'BAM',
                urls: [{'url': url}]
              }
            });
          } else {
            let url = `http://localhost/sample?accession=${id}&format=BAM`;
            let values = JSON.parse(body);
            expect(values.htsget).toBeDefined();
            expect(values.htsget.urls).toBeDefined();
            expect(values.htsget.format).toEqual("BAM");
            expect(values.htsget.urls[0].url).toEqual(url);
            expect(values.htsget.urls[0].headers).toEqual({"encoded_regions": "H4sIAAAAAAAAA4uuVipKTUstSs1LTvVLzE1VslJKzigyVNJRSs1LUbIyqY0FAMJeFOgiAAAA"});
          }
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam", "regions":[{"referenceName":"chr1", "end":4}]}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`successful redirection, range start and end given ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      if (obj.method === "GET") {
        obj.path = server_path + '?referenceName=chr1&start=4&end=400';
      } else {
        obj.path = server_path;
      }
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          expect(response.headers['content-type']).toMatch(
            /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
          );
          expect(response.statusCode).toEqual(200);
          expect(response.statusMessage).toEqual(
            'OK, see redirection instructions in the body of the message');
          if (obj.method === "GET") {
            let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1%3A5-400`;
            expect(JSON.parse(body)).toEqual({
              htsget: {
                format: 'BAM',
                urls: [{'url': url}]
              }
            });
          } else {
            let url = `http://localhost/sample?accession=${id}&format=BAM`;
            let values = JSON.parse(body);
            expect(values.htsget).toBeDefined();
            expect(values.htsget.urls).toBeDefined();
            expect(values.htsget.format).toEqual("BAM");
            expect(values.htsget.urls[0].url).toEqual(url);
            expect(values.htsget.urls[0].headers).toEqual({"encoded_regions": "H4sIAAAAAAAAA4uuVipKTUstSs1LTvVLzE1VslJKzigyVNJRKi5JLCpRsjLRUUrNSwHSBga1sQD/h2zlLgAAAA=="});
          }
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam", "regions":[{"referenceName":"chr1", "start":4, "end":400}]}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    ['bam', 'BAM', 'sam', 'SAM', 'cram', 'CRAM', 'vcf', 'VCF'].forEach( ( value ) => {
      it(`successful redirection, query with all possible params for ${optionsList.method} and ${value}`, function(done) {
        let obj = JSON.parse(JSON.stringify(optionsList));
        if (obj.method === "GET") {
          obj.path = server_path + `?referenceName=chr1&start=4&end=400&format=${value}`;
        } else {
          obj.path = server_path;
        }
        let req = http.request( obj, function(response) {
          var body = '';
          response.on('data', function(d) { body += d;});
          response.on('end', function() {
            expect(response.headers['content-type']).toMatch(
              /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
            );
            expect(response.statusCode).toBe(200);
            expect(response.statusMessage).toBe(
              'OK, see redirection instructions in the body of the message');
            let formatUpperCase = value.toUpperCase();
            if (obj.method === "GET") {
              let url = `http://localhost/sample?accession=${id}&format=${formatUpperCase}&region=chr1%3A5-400`;
              expect(JSON.parse(body)).toEqual({
                htsget: {
                  format: `${formatUpperCase}`,
                  urls: [{'url': url}]
                }
              });
            } else {
              let url = `http://localhost/sample?accession=${id}&format=${formatUpperCase}`;
              let values = JSON.parse(body);
              expect(values.htsget).toBeDefined();
              expect(values.htsget.urls).toBeDefined();
              expect(values.htsget.format).toEqual(`${formatUpperCase}`);
              expect(values.htsget.urls[0].url).toEqual(url);
              expect(values.htsget.urls[0].headers).toEqual({"encoded_regions": "H4sIAAAAAAAAA4uuVipKTUstSs1LTvVLzE1VslJKzigyVNJRKi5JLCpRsjLRUUrNSwHSBga1sQD/h2zlLgAAAA=="});
            }
            done();
          });
        });
        if (obj.method === "POST") {
          req.write(JSON.stringify({"format":value, "regions":[{"referenceName":"chr1", "start":4, "end":400}]}), () => {
            req.end();
          });
        } else {
          req.end();
        }
      });
    });

    it(`successful redirection, unknown filter ignored for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      obj.path = server_path + '?not_a_filter=1';
      let req = http.request ( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          expect(response.headers['content-type']).toMatch(
            /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
          );
          expect(response.statusCode).toEqual(200);
          expect(response.statusMessage).toEqual(
            'OK, see redirection instructions in the body of the message');
          let url = `http://localhost/sample?accession=${id}&format=BAM`;
          expect(JSON.parse(body)).toEqual({
            htsget: {
              format: 'BAM',
              urls: [{'url': url}]
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`redirection error, range is given, reference is missing for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      if (obj.method === "GET") {
        obj.path = server_path + '?start=4&end=400';
      } else {
        obj.path = server_path;
      }
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(400);
          let errorPayload = JSON.parse(body);
          expect(errorPayload.htsget).toBeDefined();
          expect(errorPayload.htsget.error).toBe(ServerHttpError.INVALID_INPUT);
          expect(response.statusMessage).toBe(
            "'referenceName' attribute required if 'start' or 'end' attribute is given");
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam", "regions":[{"start":4, "end":400}]}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`redirection error, range start is not an integer ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      if (obj.method === "GET") {
        obj.path = server_path + '?referenceName=chr1&start=5.5&end=400';
      } else {
        obj.path = server_path;
      }
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(400);
          let errorPayload = JSON.parse(body);
          expect(errorPayload.htsget).toBeDefined();
          expect(errorPayload.htsget.error).toBe(ServerHttpError.INVALID_INPUT);
          expect(response.statusMessage).toEqual(
            "ranges must be integers");
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam", "regions":[{"referenceName":"chr1", "start":5.5, "end":400}]}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`redirection error, range start is a negative integer for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      if (obj.method === "GET") {
        obj.path = server_path + '?referenceName=chr1&start=-44&end=400';
      } else {
        obj.path = server_path;
      }
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(400);
          let errorPayload = JSON.parse(body);
          expect(errorPayload.htsget).toBeDefined();
          expect(errorPayload.htsget.error).toBe(ServerHttpError.INVALID_INPUT);
          expect(response.statusMessage).toEqual("ranges must be unsigned integers");
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam", "regions":[{"referenceName":"chr1", "start":-44, "end":400}]}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`redirection error, range end is not an integer for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      if (obj.method === "GET") {
        obj.path = server_path + '?referenceName=chr1&start=4&end=foo';
      } else {
        obj.path = server_path;
      }
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(400);
          let errorPayload = JSON.parse(body);
          expect(errorPayload.htsget).toBeDefined();
          expect(errorPayload.htsget.error).toBe(ServerHttpError.INVALID_INPUT);
          expect(response.statusMessage).toEqual("ranges must be integers");
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam", "regions":[{"referenceName":"chr1", "start":4, "end":"foo"}]}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`redirection error, range end is a negative integer for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      if (obj.method === "GET") {
        obj.path = server_path + '?referenceName=chr1&start=4&end=-400';
      } else {
        obj.path = server_path;
      }
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(400);
          let errorPayload = JSON.parse(body);
          expect(errorPayload.htsget).toBeDefined();
          expect(errorPayload.htsget.error).toBe(ServerHttpError.INVALID_INPUT);
          expect(response.statusMessage).toEqual("ranges must be unsigned integers");
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam", "regions":[{"referenceName":"chr1", "start":4, "end":-400}]}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });
    
    it(`redirection error, range start is bigger than range end for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      if (obj.method === "GET") {
        obj.path = server_path + '?referenceName=chr1&start=400&end=4';
      } else {
        obj.path = server_path;
      }
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(400);
          let errorPayload = JSON.parse(body);
          expect(errorPayload.htsget).toBeDefined();
          expect(errorPayload.htsget.error).toBe(ServerHttpError.INVALID_INPUT);
          expect(response.statusMessage).toEqual(
            'Range end should be bigger than start');
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam", "regions":[{"referenceName":"chr1", "start":400, "end":4}]}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`redirection error, unknown format requested for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      if (obj.method === "GET") {
        obj.path = server_path + '?format=fa';
      } else {
        obj.path = server_path;
      }
      let req = http.request( obj, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          checkGenericErrorResponse(expect, response);
          expect(response.statusCode).toEqual(400);
          let errorPayload = JSON.parse(body);
          expect(errorPayload.htsget).toBeDefined();
          expect(errorPayload.htsget.error).toEqual(ServerHttpError.UNSUPPORTED_FORMAT);
          expect(response.statusMessage).toEqual(
            "Format 'fa' is not supported, supported formats: BAM, CRAM, SAM, VCF");
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"fa"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });
  });
});

describe('redirection when running behind a proxy', () => {
  const server = http.createServer();
  let socket = tmp.tmpNameSync();
  let id              = 'EGA45678';
  let serverPath      = `${GA4GH_URL}/${id}`;

  beforeAll((done) =>  {
    fse.ensureDirSync(tmpDir);
    config.provide(() => {return {
      tempdir:   tmpDir,
      skipauth:  true,
      proxylist: {
        'http://myserver.com':      'http://myserver.com/path1/path2',
        'http://myserver.com:3456': 'http://myserver.com:3456/path3',
      }
    };});
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {});
      c.handleRequest();
    });
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });

  afterAll( () => {
    server.close();
    utils.removeSocket(socket);
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });

  let requestOptions = [
    {method: 'GET', socketPath: socket, headers: {}},
    {method: 'POST', socketPath: socket, headers: {}}
  ];

  requestOptions.forEach( optionsList => {
    it(`direct access is not allowed - GA4GH url for ${optionsList.method}`, (done) => {
      let obj = JSON.parse(JSON.stringify(optionsList));
      obj.path = '/sample/' + id;
      let req = http.request( obj );
      req.on('response', (res) => {
        expect(res.statusCode).toEqual(403);
        expect(res.statusMessage).toEqual(
          'Bypassing proxy server is not allowed');
        done();
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`direct access is not allowed - sample url for ${optionsList.method}`, (done) => {
      let obj = JSON.parse(JSON.stringify(optionsList));
      obj.path = serverPath;
      let req = http.request( obj );
      req.on('response', (res) => {
        expect(res.statusCode).toEqual(403);
        expect(res.statusMessage).toEqual(
          'Bypassing proxy server is not allowed');
        done();
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`unknown proxy is not allowed for ${optionsList.method}`, (done) => {
      let obj = JSON.parse(JSON.stringify(optionsList));
      obj.path = serverPath;
      obj.headers =  {'X-Forwarded-Host': 'myserver.com:9090'};
      let req = http.request( obj );
      req.on('response', (res) => {
        expect(res.statusCode).toEqual(403);
        expect(res.statusMessage).toEqual(
          'Unknown proxy http://myserver.com:9090');
        done();
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`Redirection to one of known proxies for ${optionsList.method}`, (done) => {
      let obj = JSON.parse(JSON.stringify(optionsList));
      obj.path = serverPath;
      obj.headers =  {'X-Forwarded-Host': 'myserver.com:3456'};
      let url = `http://myserver.com:3456/path3/sample?accession=${id}&format=BAM`;
      let req = http.request( obj );
      req.on('response', (res) => {
        var body = '';
        expect(res.statusCode).toEqual(200);
        res.on('data', (d) => { body += d;});
        res.on('end', () => {
          expect(JSON.parse(body)).toEqual({
            htsget: {
              format: 'BAM',
              urls: [{'url': url}]
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    it(`Redirection to one of known proxies for ${optionsList.method}`, (done) => {
      let obj = JSON.parse(JSON.stringify(optionsList));
      obj.path = serverPath;
      obj.headers =  {'X-Forwarded-Host': 'myserver.com'};
      let url = `http://myserver.com/path1/path2/sample?accession=${id}&format=BAM`;
      let req = http.request( obj );
      req.on('response', (res) => {
        var body = '';
        expect(res.statusCode).toEqual(200);
        res.on('data', (d) => { body += d;});
        res.on('end', () => {
          expect(JSON.parse(body)).toEqual({
            htsget: {
              format: 'BAM', urls: [{'url': url}]
            }
          });
          done();
        });
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });   
  });
});

describe('content type', function() {
  const server = http.createServer();
  var socket = tmp.tmpNameSync();

  beforeAll((done) => {
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });
  afterAll( () => {
    server.close();
    utils.removeSocket(socket);
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });

  let requestOptions = [
    {method: 'GET', socketPath: socket, path: '/file'},
    {method: 'POST', socketPath: socket, path: '/file'}
  ];
  
  requestOptions.forEach( optionsList => {
    it(`data format driven content type for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      server.on('request', (request, response) => {
        let c = new RangerController(request, response, {one: "two"});
        expect( () => {c.contentType();} )
          .toThrowError(assert.AssertionError,
                        'Non-empty format string should be given');
        expect(c.contentType('SAM')).toBe('text/vnd.ga4gh.sam');
        expect(c.contentType('VCF')).toBe('text/vnd.ga4gh.vcf');
        expect(c.contentType('BAM')).toBe('application/vnd.ga4gh.bam');
        expect(c.contentType('CRAM')).toBe('application/vnd.ga4gh.cram');
        done();
      });

      let req = http.request( obj, function(response) {
        let body = '';
        response.on('data', function(d) { body += d;});
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });
  });
});

describe('trailers in response', function() {
  const server = http.createServer();
  var socket = tmp.tmpNameSync();

  beforeAll((done) => {
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });
  afterAll( () => {
    server.close();
    utils.removeSocket(socket);
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });
  
  let requestOptions = [
    {method: 'GET', socketPath: socket, path: '/file', headers: {}},
    {method: 'POST', socketPath: socket, path: '/file', headers: {}}
  ];

  requestOptions.forEach ( optionsList => {
    it(`no trailers without TE header for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        let c = new RangerController(request, response, {one: "two"});
        expect(c.sendTrailer).toBe(false);
        done();
      });
      let req = http.request( obj, function() {});
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });

    ['TE', 'te', 'Te', 'tE'].forEach( ( headerName ) => {
      it(`trailers with ${headerName} header for ${optionsList.method}`, function(done) {
        let obj = JSON.parse(JSON.stringify(optionsList));
        server.removeAllListeners('request');
        server.on('request', (request, response) => {
          let c = new RangerController(request, response, {one: "two"});
          expect(c.sendTrailer).toBe(true);
          done();
        });
        let headers = {};
        headers[headerName] = 'trailers';
        obj.headers = headers;
        let req = http.request( obj, () => {});
        if (obj.method === "POST") {
          req.write(JSON.stringify({"format":"bam"}), () => {
            req.end();
          });
        } else {
          req.end();
        }
      });
    });
  });
});

describe('CORS in response', function() {
  var server;
  let serverPath = `${GA4GH_URL}/EGA45678`;
  let socket = tmp.tmpNameSync();

  let checkHeaders = (headers, origin) => {
    expect(headers.vary).toBe('Origin', 'Vary header is set');
    expect(headers['access-control-allow-origin']).toBe(
      origin, `allowed origin is ${origin}`);
    expect(headers['access-control-allow-methods']).toBe(
      'GET,OPTIONS', 'allowed methods are set');
    expect(headers['access-control-allow-headers']).toBe(
      'TE,X-Remote-User,withcredentials', 'allowed headers are set');
    expect(headers['access-control-max-age']).toBe('1800', 'max age is set');
    expect(Object.keys(headers).indexOf('Access-Control-Allow-Credentials')).toBe(
      -1, 'Access-Control-Allow-Credentials header is not set');
  };

  beforeAll((done) => {
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
    server = http.createServer();
    server.listen(socket, () => { console.log('listening'); done();});
  });
  afterAll( () => {
    server.close();
    utils.removeSocket(socket);
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });

  let requestOptions = [
    {method: 'GET', socketPath: socket, path: serverPath},
    {method: 'POST', socketPath: socket, path: serverPath}
  ];
  
  requestOptions.forEach( optionsList => {
    it(`no CORS in a response to a standart request for ${optionsList.method}`, function(done) {
      let obj = JSON.parse(JSON.stringify(optionsList));
      config.provide( () => {
        return {tempdir: tmpDir, anyorigin: false, originlist: null, skipauth: true};
      });
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        expect('origin' in request.headers).toBe(false, 'request does not have Origin header');
        let c = new RangerController(request, response, {one: "two"});
        c.handleRequest();
      });

      let req = http.request( obj, (res) => {
        expect( Object.keys(res.headers).filter((headerName) => {
          return headerName.startsWith('Access-Control');
        }).length).toBe(0, 'no CORS headers in reply');
        expect(res.headers.vary).toBe('Origin', 'Vary header is set');
        done();
      });
      if (obj.method === "POST") {
        req.write(JSON.stringify({"format":"bam"}), () => {
          req.end();
        });
      } else {
        req.end();
      }
    });
  });

  it('no CORS headers in a response to CORS GET request due to server options', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: false, originlist: null, skipauth: true};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect('origin' in request.headers).toBe(true, 'request has Origin header');
      let c = new RangerController(request, response, {one: "two"});
      c.handleRequest();
    });

    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'GET'};
    let req = http.request(options);
    req.on('response', (res) => {
      expect( Object.keys(res.headers).filter((headerName) => {
        return headerName.startsWith('access-control');
      }).length).toBe(0, 'no CORS headers in reply');
      expect(res.headers.vary).toBe('Origin', 'Vary header is set');
      done();
    });
    req.end();
  });

  it('no CORS headers in a response to CORS OPTIONS request due to server options', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: false, originlist: null, skipauth: true};
    });
    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'OPTIONS'};
    let req = http.request(options);
    req.on('response', (res) => {
      expect( Object.keys(res.headers).filter((headerName) => {
        return headerName.startsWith('access-control');
      }).length).toBe(0, 'no CORS headers in reply');
      done();
    });
    req.end();
  });
  
  it('no CORS headers in a response to CORS POST request due to server options', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: false, originlist: null, skipauth: true};
    });
    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'POST'};
    let req = http.request(options);
    req.on('response', (res) => {
      expect( Object.keys(res.headers).filter((headerName) => {
        return headerName.startsWith('access-control');
      }).length).toBe(0, 'no CORS headers in reply');
      done();
    });
    req.write(JSON.stringify({"format":"bam"}), () => {
      req.end();
    });
  });

  it('Allow all CORS headers in a response to CORS GET request', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: true, originlist: null, skipauth: true};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect('origin' in request.headers).toBe(true, 'request has Origin header');
      let c = new RangerController(request, response, {one: "two"});
      c.handleRequest();
    });

    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'GET'};
    let req = http.request(options);
    req.on('response', (res) => {
      checkHeaders(res.headers, '*');
      done();
    });
    req.end();
  });

  it('Allow all CORS headers in a response to CORS OPTIONS request', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: true, originlist: null, skipauth: true};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect('origin' in request.headers).toBe(true, 'request has Origin header');
      let c = new RangerController(request, response, {one: "two"});
      c.handleRequest();
    });

    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'OPTIONS'};
    let req = http.request(options);
    req.on('response', (res) => {
      checkHeaders(res.headers, '*');
      done();
    });
    req.end();
  });
  
  it('Allow all CORS headers in a response to CORS POST request', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: true, originlist: null, skipauth: true};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect('origin' in request.headers).toBe(true, 'request has Origin header');
      let c = new RangerController(request, response, {one: "two"});
      c.handleRequest();
    });

    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'POST'};
    let req = http.request(options);
    req.on('response', (res) => {
      checkHeaders(res.headers, '*');
      done();
    });
    req.write(JSON.stringify({"format":"bam"}), () => {
      req.end();
    });
  });

  it('No CORS headers since the origin is not white listed', function(done) {
    config.provide( () => {
      return {tempdir:     tmpDir,
              anyorigin:   false,
              originlist: ['http://other.com','http://other.com:9090'],
              skipauth:   true};
    });
    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'OPTIONS'};
    let req = http.request(options);
    req.end();
    req.on('response', (res) => {
      expect( Object.keys(res.headers).filter((headerName) => {
        return headerName.startsWith('access-control');
      }).length).toBe(0, 'no CORS headers in reply');
      done();
    });
  });

  it('Origin-specific CORS headers set for a white listed origin', function(done) {
    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://other.com'},
                   method:     'OPTIONS'};
    let req = http.request(options);
    req.on('response', (res) => {
      checkHeaders(res.headers, 'http://other.com');
      done();
    });
    req.end();
  });

  it('Additional CORS header is set when running with authorization', function(done) {
    config.provide( () => {
      return {tempdir:    tmpDir,
              anyorigin:  false,
              originlist: ['http://other.com','http://other.com:9090'],
              skipauth:   false};
    });
    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://other.com:9090'},
                   method:    'OPTIONS'};
    let req = http.request(options);
    req.on('response', (res) => {
      expect(res.headers['access-control-allow-origin']).toBe(
        'http://other.com:9090', 'origin http://other.com:9090 is allowed');
      expect(res.headers['access-control-allow-credentials']).toBe(
        'true', 'Access-Control-Allow-Credentials header is set');
      done();
    });
    req.end();
  });
});

describe('POST regions in header', () => {
  const server = http.createServer();
  var socket   = tmp.tmpNameSync();
  let id       = 'EGA45678';
  let server_path = `${GA4GH_TOKEN_URL}/${id}`;

  beforeAll((done) =>  {
    fse.ensureDirSync(tmpDir);
    options = config.provide(() => {
      return { tempdir: tmpDir };
    });
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {});
      c.handleRequest();
    });
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });

  afterAll( () => {
    server.close();
    utils.removeSocket(socket);
    try { fse.removeSync(tmpDir); } catch (e) { console.log(e); }
  });

  
  let obj =  {method: 'POST', socketPath: socket, path: server_path};
  
  it('Three regions non-overlap', function(done) {
    let my_token = 'XXXYYYXXX';
    let headers  = {};
    headers[
      constants.TOKEN_BEARER_KEY_NAME
    ] = tokenUtils.formatTokenForHeader(my_token);
    obj.headers = headers;
    let req = http.request( obj, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toMatch(
          /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
        );
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM`;
        expect(JSON.parse(body)).toEqual({
          htsget: {
            format: 'BAM',
            urls: [
              {
                'url': url,
                'headers': {
                  'Authorization': `Bearer ${my_token}`,
                  'encoded_regions': 'H4sIAAAAAAAAA4uuVipKTUstSs1LTvVLzE1VslJKzigyVKrVwSphpKSjVFySWFSiZGVpaamjlJqXomRlaGBgQFi9EVAVVIMRUEdtLAAiqs+mewAAAA=='
                  // [{"referenceName":"chr1"},{"referenceName":"chr2","start":999,"end":1000},{"referenceName":"chr2","start":2000,"end":2100}]
                }
              }
            ]
          }
        });
        done();
      });
    });
    if (obj.method === "POST") {
      req.write(JSON.stringify({"format":"bam","regions" : [
        { "referenceName" : "chr1" },
        { "referenceName" : "chr2", "start" : 999, "end" : 1000 },
        { "referenceName" : "chr2", "start" : 2000, "end" : 2100 }
      ] }), () => {
        req.end();
      });
    } else {
      req.end();
    }
  });

  obj =  {method: 'POST', socketPath: socket, path: server_path};
  
  it('Five regions with overlap merge', function(done) {
    let my_token = 'XXXYYYXXX';
    let headers  = {};
    headers[
      constants.TOKEN_BEARER_KEY_NAME
    ] = tokenUtils.formatTokenForHeader(my_token);
    obj.headers = headers;
    let req = http.request( obj, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toMatch(
          /application\/vnd\.ga4gh\.htsget\.\S+\+json/i
        );
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM`;
        expect(JSON.parse(body)).toEqual({
          htsget: {
            format: 'BAM',
            urls: [
              {
                'url': url,
                'headers': {
                  'Authorization': `Bearer ${my_token}`,
                  'encoded_regions': 'H4sIAAAAAAAAA4uuVipKTUstSs1LTvVLzE1VslJKzigyVKrVwSphpKSjVFySWFSiZGVkYqqjlJqXomRlZm5hSli9oZGBiQlUh5GhgYEZDi3GCC2mZsYGtbEAs5Ckf6MAAAA='
                  // [{"referenceName":"chr1"},{"referenceName":"chr2","start":245,"end":6785},{"referenceName":"chr2","start":12044,"end":21006},{"referenceName":"chr3","start":5630}]
                }
              }
            ]
          }
        });
        done();
      });
    });
    if (obj.method === "POST") {
      req.write(JSON.stringify({"format":"bam","regions" : [
        { "referenceName" : "chr1" },
        { "referenceName" : "chr2", "start" : 245, "end" : 6785 },
        { "referenceName" : "chr2", "start" : 12044, "end" : 21006 },
        { "referenceName" : "chr3", "start" : 9200},
        { "referenceName" : "chr3", "start" : 5630, "end" : 12312}
      ] }), () => {
        req.end();
      });
    } else {
      req.end();
    }
  });
  
});