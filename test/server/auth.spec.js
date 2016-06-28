/* globals describe, it, expect, beforeAll, afterAll */

"use strict";

const assert      = require('assert');
const child       = require('child_process');
const MongoClient = require('mongodb').MongoClient;
const tmp         = require('tmp');
const fse         = require('fs-extra');

const DataAccess  = require('../../lib/server/auth.js');

const BASE_PORT  = 1100;
const PORT_RANGE = 200;
const PORT = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;
const FIXTURES = 'test/server/data/fixtures/access_control_group.json';

describe('Authorisation', function() {
  var tmpobj = tmp.dirSync({ prefix: 'npg_ranger_test_' });
  var tmp_dir = tmpobj.name;
  console.log(`MONGO data directory: ${tmp_dir}`);
  var db_name = 'npg_ranger_test';
  var url = `mongodb://localhost:${PORT}/${db_name}`;

  beforeAll( () => {
    let command = `mongod -f test/server/data/mongodb_conf.yml --port ${PORT} --dbpath ${tmp_dir} --pidfilepath ${tmp_dir}/mpid --logpath ${tmp_dir}/dbserver.log`;
    console.log(`\nCommand to start MONGO DB daemon: ${command}`);
    let out = child.execSync(command);
    console.log(`Started MONGO DB daemon: ${out}`);
    command = `mongoimport --port ${PORT} --db ${db_name} --collection access_control_group --jsonArray --file ${FIXTURES}`;
    out = child.execSync(command);
    console.log(`Loaded data to MONGO DB: ${out}`);
  });

  it('Input validation', function() {
    expect( () => {new DataAccess();} ).toThrowError(ReferenceError,
      'Database handle is required');
    expect( () => {
      let da = new DataAccess({});
      da.authorise();
    }).toThrowError(ReferenceError, 'Username is required');
    expect( () => {
      let da=new DataAccess({});
      da.authorise('alice');
    }).toThrowError(Error, 'Access groups array is not available');
    expect( () => {
      let da=new DataAccess({});
      da.authorise('alice', 2);
    }).toThrowError(Error, 'Access groups array is not available');
    expect( () => {
      let da=new DataAccess({});
      da.authorise('alice', []);
    }).toThrowError(Error, 'Access groups array is not available');
  });

  it('Authorisation failed - access_control_group_id value is missing', function(done) {
    var da = new DataAccess({});
    da.on('failed', (username, reason) => {
      expect(username).toBe('alice');
      expect(reason).toBe(
        'Some access group ids are not defined');
      done();
    });
    da.authorise('alice', ["9", "", "10"]);
  });

  it('Authorised - files belong to the same auth group', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var da = new DataAccess(db);
      da.on('authorised', (username) => {
         expect(username).toBe('alice');
         done();
      });
      // We are not interested in what happens on fail.
      // If 'authorised' event is not processed withing the set
      // time limit, the test will fail.
      da.authorise('alice', ["6", "6"]);
    });
  });

  it('Authorised - files belong to different auth group', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var da = new DataAccess(db);
      da.on('authorised', (username) => {
         expect(username).toBe('alice');
         done();
      });
      da.authorise('alice', ["7", "6"]);
    });
  });

  it('Authorisation failed - no auth for some of the files', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var da = new DataAccess(db);
      da.on('failed', (username, reason) => {
         expect(username).toBe('alice');
         expect(reason).toBe('Not authorised for some of the files');
         done();
      });
      da.authorise('alice', ["8", "6"]);
    });
  });

  it('Authorisation failed - no auth for any of the files', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var da = new DataAccess(db);
      da.on('failed', (username, reason) => {
         expect(username).toBe('alice');
         expect(reason).toBe('Not authorised for any of the files');
         done();
      });
      da.authorise('alice', ["8", "9"]);
    });
  });

  afterAll(() => {
    let out = child.execSync(`mongod --shutdown --dbpath ${tmp_dir} --pidfilepath ${tmp_dir}/mpid`);
    console.log(`\nMONGODB server has been shut down: ${out}`);
    fse.remove(tmp_dir, (err) => {if (err) {console.log(`Error removing ${tmp_dir}: ${err}`);}});
  });
});

