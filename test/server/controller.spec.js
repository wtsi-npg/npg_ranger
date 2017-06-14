/* globals describe, it, expect, beforeAll, afterAll, jasmine */

"use strict";

const assert  = require('assert');
const http    = require('http');
const fs      = require('fs');
const fse     = require('fs-extra');
const tmp     = require('tmp');
const RangerController = require('../../lib/server/controller.js');
const config  = require('../../lib/config.js');

const utils   = require('./test_utils.js');

// Create temp dir here so it is available for all tests.
// Use this dir as a default dir that will be available in all.
var tmpDir    = config.tempFilePath('npg_ranger_controller_test_');
var dummy     = function() { return {tempdir: tmpDir, skipauth: true}; };
var options;

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

  it('db object is not given or is not an object - error', function(done) {

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
      response.end();
      done();
    });
    http.get({socketPath: socket}, function() {});
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

  it('Method not allowed error', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });

    let options = { socketPath: socket, method: 'POST' };
    let req = http.request(options);
    req.on('response', (response) => {
      expect(response.headers['content-type']).toEqual('application/json');
      expect(response.statusCode).toEqual(405);
      expect(response.statusMessage).toEqual('POST request is not allowed');
      done();
    });
    req.end();
  });

  it('Authentication error', function(done) {
    config.provide( () => {
      return {tempdir:  tmpDir,
              skipauth: false};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });
    let httpopts = {
      socketPath: socket,
      path: '/authuser',
    };
    http.get(httpopts, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(401);
        expect(response.statusMessage).toEqual('Proxy authentication required');
        expect(JSON.parse(body)).toEqual({
          error:   "InvalidAuthentication",
          message: "Proxy authentication required"
         });
        done();
      });
    });
  });

  it('Not found error, no auth', function(done) {
    config.provide( () => {
      return {tempdir:  tmpDir,
              skipauth: true};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });
    http.get({socketPath: socket, path: '/invalid'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : /invalid');
        expect(JSON.parse(body)).toEqual({
          error:   "NotFound",
          message: "URL not found : /invalid"
        });
        done();
      });
    });
  });

  it('Not found error, auth checked', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });
    let req = http.request({socketPath: socket, path: '/invalid'});
    req.setHeader('X-Remote-User', 'user1');
    req.on('response', function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : /invalid');
        expect(JSON.parse(body)).toEqual({
          error:   "NotFound",
          message: "URL not found : /invalid"
        });
        done();
      });
    });
    req.end();
  });

  it('Invalid input error for duplicate values for an attribute', ( done ) => {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });

    http.get({socketPath: socket, path: '/sample?attr1=value1&attr2=value2&attr1=value3'},
      ( response ) => {
      var body = '';
      response.on('data', ( d ) => { body += d;});
      response.on('end', () => {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        let m = "Invalid request: multiple values for attribute 'attr1'";
        expect(response.statusMessage).toEqual(m);
        expect(JSON.parse(body)).toEqual({
          error:   "InvalidInput",
          message: m
        });
        done();
      });
    });
  });


  it('Invalid input error for a sample url', ( done ) => {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });

    http.get({socketPath: socket, path: '/sample'}, ( response ) => {
      var body = '';
      response.on('data', ( d ) => { body += d;});
      response.on('end', () => {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        let m = 'Invalid request: sample accession number should be given';
        expect(response.statusMessage).toEqual(m);
        expect(JSON.parse(body)).toEqual({
          error:   "InvalidInput",
          message: m
        });
        done();
      });
    });
  });

  it('Invalid input error for a file url', ( done ) => {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });

    http.get({socketPath: socket, path: '/file'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        let m = 'Invalid request: file name should be given';
        expect(response.statusMessage).toEqual(m);
        expect(JSON.parse(body)).toEqual({
          error:   "InvalidInput",
          message: m
        });
        done();
      });
    });
  });

  it('Invalid input error for a vcf file when multiref set', (done) => {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      // Set multiref mode
      config.provide(() => { return {tempdir: tmpDir, multiref: true, skipauth: true}; });
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
      // unset multiref
      config.provide(dummy);
    });

    http.get({socketPath: socket, path: '/sample?accession=XYZ120923&format=vcf'}, (response) => {
      var body = '';
      response.on('data', (d) => { body += d;});
      response.on('end', () => {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        let m = 'Invalid request: cannot produce VCF files while multiref set on server';
        expect(response.statusMessage).toEqual(m);
        expect(JSON.parse(body)).toEqual({
          error:   "InvalidInput",
          message: m
        });
        done();
      });
    });
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
          expect(response.headers['content-type']).toEqual('application/json');
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
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual(
          `No files for sample accession ${missingAcc}`
        );
        done();
      });
    });
  });
});

describe('Redirection in json response', function() {
  const server = http.createServer();
  var socket   = tmp.tmpNameSync();
  let id       = 'EGA45678';
  let server_path_basic = '/ga4gh/v.0.1/get/sample';
  let server_path = server_path_basic + '/' + id;

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

  it('invalid url - no id - error response', function(done) {
    http.get({socketPath: socket, path: server_path_basic}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : ' + server_path_basic);
        done();
      });
    });
  });

  it('invalid url - no id - error response', function(done) {
    let path = server_path_basic + '/';
    http.get({socketPath: socket, path: path}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : ' + path);
        done();
      });
    });
  });

  it('invalid sample id - error response', function(done) {
    let path = server_path_basic + 'ERS-4556';
    http.get({socketPath: socket, path: path}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : ' + path);
        done();
      });
    });
  });

  it('successful redirection, no query params', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('successful redirection, format given', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?format=CRAM'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=CRAM`;
        expect(JSON.parse(body)).toEqual({format: 'CRAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('successful redirection, chromosome given', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('successful redirection, range start given', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=3'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1%3A4`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('successful redirection, range end given', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&end=4'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1%3A1-4`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('successful redirection, range start and end given', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=4&end=400'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1%3A5-400`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  ['bam', 'BAM', 'sam', 'SAM', 'cram', 'CRAM', 'vcf', 'VCF'].forEach( ( value ) => {
    it('successful redirection, query with all possible params', function(done) {
      http.get(
        { socketPath: socket,
          path: server_path + `?referenceName=chr1&start=4&end=400&format=${value}`}, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          expect(response.headers['content-type']).toEqual('application/json');
          expect(response.statusCode).toBe(200);
          expect(response.statusMessage).toBe(
            'OK, see redirection instructions in the body of the message');
          let formatUpperCase = value.toUpperCase();
          let url = `http://localhost/sample?accession=${id}&format=${formatUpperCase}&region=chr1%3A5-400`;
          expect(JSON.parse(body)).toEqual({format: `${formatUpperCase}`, urls: [{'url': url}]});
          done();
        });
      });
    });
  });

  it('successful redirection, filter given', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?target=0&manual_qc=&alignment_not=undef'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM&target=0&manual_qc=&alignment_not=undef`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('successful redirection, unknown filter ignored', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?not_a_filter=1'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('redirection error, range is given, reference is missing', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?start=4&end=400'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        expect(response.statusMessage).toBe(
          "'referenceName' attribute requered if 'start' or 'end' attribute is given");
        done();
      });
    });
  });

  it('redirection error, range start is not an integer', function(done) {
    http.get(
     {socketPath: socket,
      path: server_path + '?referenceName=chr1&start=5.5&end=400'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        expect(response.statusMessage).toEqual(
          "'5.5' is not an integer");
        done();
      });
    });
  });

  it('redirection error, range start is a negative integer', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=-44&end=400'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        expect(response.statusMessage).toEqual("'-44' is not an unsigned integer");
        done();
      });
    });
  });

  it('redirection error, range end is not an integer', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=4&end=foo'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        expect(response.statusMessage).toEqual("'foo' is not an integer");
        done();
      });
    });
  });

  it('redirection error, range end is a negative integer', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=4&end=-400'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        expect(response.statusMessage).toEqual("'-400' is not an unsigned integer");
        done();
      });
    });
  });

  it('redirection error, range start is bigger than range end', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=400&end=4'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        expect(response.statusMessage).toEqual(
          'Range end should be bigger than start');
        done();
      });
    });
  });

  it('redirection error, invalid characted in reference name', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr@1&start=400&end=4'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        expect(response.statusMessage).toEqual(
          'Invalid character in reference name chr@1');
        done();
      });
    });
  });

  it('redirection error, unknown format requested', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?format=fa'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(400);
        expect(response.statusMessage).toEqual(
          "Format 'fa' is not supported, supported formats: BAM, CRAM, SAM, VCF");
        done();
      });
    });
  });

});

describe('redirection when running behind a proxy', () => {
  const server = http.createServer();
  let socket = tmp.tmpNameSync();
  let id              = 'EGA45678';
  let serverPath      = '/ga4gh/v.0.1/get/sample/' + id;

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

  it('direct access is not allowed - GA4GH url', (done) => {
    let options = {
      socketPath: socket,
      path:       '/sample/' + id,
      headers:    {},
      method:    'GET'};
    let req = http.request(options);
    req.on('response', (res) => {
      expect(res.statusCode).toEqual(403);
      expect(res.statusMessage).toEqual(
        'Bypassing proxy server is not allowed');
      done();
    });
    req.end();
  });

  it('direct access is not allowed - sample url', (done) => {
    let options = {
      socketPath: socket,
      path:       serverPath,
      headers:    {},
      method:    'GET'};
    let req = http.request(options);
    req.on('response', (res) => {
      expect(res.statusCode).toEqual(403);
      expect(res.statusMessage).toEqual(
        'Bypassing proxy server is not allowed');
      done();
    });
    req.end();
  });

  it('unknown proxy is not allowed', (done) => {
    let options = {
      socketPath: socket,
      path:       serverPath,
      headers:    {'X-Forwarded-Host': 'myserver.com:9090'},
      method:    'GET'};
    let req = http.request(options);
    req.on('response', (res) => {
      expect(res.statusCode).toEqual(403);
      expect(res.statusMessage).toEqual(
        'Unknown proxy http://myserver.com:9090');
      done();
    });
    req.end();
  });

  it('Redirection to one of known proxies', (done) => {
    let url = `http://myserver.com:3456/path3/sample?accession=${id}&format=BAM`;
    let options = {
      socketPath: socket,
      path:       serverPath,
      headers:    {'X-Forwarded-Host': 'myserver.com:3456'},
      method:    'GET'};
    let req = http.request(options);
    req.on('response', (res) => {
      var body = '';
      expect(res.statusCode).toEqual(200);
      res.on('data', (d) => { body += d;});
      res.on('end', () => {
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
    req.end();
  });

  it('Redirection to one of known proxies', (done) => {
    let url = `http://myserver.com/path1/path2/sample?accession=${id}&format=BAM`;
    let options = {
      socketPath: socket,
      path:       serverPath,
      headers:    {'X-Forwarded-Host': 'myserver.com'},
      method:    'GET'};
    let req = http.request(options);
    req.on('response', (res) => {
      var body = '';
      expect(res.statusCode).toEqual(200);
      res.on('data', (d) => { body += d;});
      res.on('end', () => {
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
    req.end();
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

  it('data format driven content type', function(done) {
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

    http.get({socketPath: socket, path: '/file'}, function(response) {
      let body = '';
      response.on('data', function(d) { body += d;});
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

  it('no trailers without TE header', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect(c.sendTrailer).toBe(false);
      done();
    });

    http.get({socketPath:socket, path: '/file'}, function() {});
  });

  ['TE', 'te', 'Te', 'tE'].forEach( ( headerName ) => {
    it(`trailers with ${headerName} header`, function(done) {
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        let c = new RangerController(request, response, {one: "two"});
        expect(c.sendTrailer).toBe(true);
        done();
      });

      let headers = {};
      headers[headerName] = 'trailers';
      http.get({socketPath:socket, path: '/file', headers: headers}, () => {});
    });
  });
});

describe('CORS in response', function() {
  var server;
  let serverPath = '/ga4gh/v.0.1/get/sample/EGA45678';
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

  it('no CORS in a response to a standart request', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: false, originlist: null, skipauth: true};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect('origin' in request.headers).toBe(false, 'request does not have Origin header');
      let c = new RangerController(request, response, {one: "two"});
      c.handleRequest();
    });

    http.get({socketPath:socket, path: serverPath}, (res) => {
      expect( Object.keys(res.headers).filter((headerName) => {
        return headerName.startsWith('Access-Control');
      }).length).toBe(0, 'no CORS headers in reply');
      expect(res.headers.vary).toBe('Origin', 'Vary header is set');
      done();
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
