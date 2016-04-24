"use strict";

const assert      = require('assert');
const child       = require('child_process');
const MongoClient = require('mongodb').MongoClient;
const tmp         = require('tmp');
const fse         = require('fs-extra');

const DataAccess  = require('../lib/auth.js');

jasmine.DEFAULT_TIMEOUT_INTERVAL = 5000;

const BASE_PORT  = 1100;
const PORT_RANGE = 200;
const PORT = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;
const FIXTURES = 'test/data/fixtures/access_control_group.json';

describe('Authorisation', function() {
  var tmpobj = tmp.dirSync({ prefix: 'npg_ranger_test_' });
  var tmp_dir = tmpobj.name;
  console.log(`MONGO data directory: ${tmp_dir}`);
  var db_name = 'npg_ranger_test';
  var url = `mongodb://localhost:${PORT}/${db_name}`;

  beforeAll( () => {
    let command = `mongod -f test/data/mongodb_conf.yml --port ${PORT} --dbpath ${tmp_dir} --pidfilepath ${tmp_dir}/mpid --logpath ${tmp_dir}/dbserver.log`;
    console.log(`\nCommand to start MONGO DB daemon: ${command}`);
    let out = child.execSync(command);
    console.log(`Started MONGO DB daemon: ${out}`);
    command = `mongoimport --port ${PORT} --db ${db_name} --collection access_control_group --jsonArray --file ${FIXTURES}`;
    out = child.execSync(command);
    console.log(`Loaded data to MONGO DB: ${out}`);
  });

  it('Authorised - files belong to the same auth group', function(done) {
    
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var files = [
         {"filepath_by_host" : {
     		"irods-seq-i12" : "/irods-seq-i12-de/seq/foo",
     		"irods-seq-sr02" : "/irods-seq-sr02-ddn-ra08-0-1-2/seq/foo",
     		"*" : "irods:/seq/foo1"
     	                      }, "access_control_group_id": "6"},
     	{"filepath_by_host" : {
     		"irods-seq-i12" : "/irods-seq-i12-de/seq/foo",
     		"irods-seq-sr02" : "/irods-seq-sr02-ddn-ra08-0-1-2/seq/foo",
     		"*" : "irods:/seq/foo2"
     	                      }, "access_control_group_id": "6"}
                  ];
      
      var da = new DataAccess(db, files);
      da.on('authorised', (username) => {
         expect(username).toBe('alice');
         done();
      });
      // We are not interested in what happens on fail.
      // If 'authorised' event is not processed withing the set
      // time limit, the test will fail.
      da.authorise('alice');
    });
  });

  it('Authorised - files belong to different auth group', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var files = [
        {"filepath_by_host" : {"*" : "irods:/seq/foo1"},
          "access_control_group_id": "6"},
        {"filepath_by_host" : {"*" : "irods:/seq/foo2"},
          "access_control_group_id": "7"}
                  ];
      var da = new DataAccess(db, files);
      da.on('authorised', (username) => {
         expect(username).toBe('alice');
         done();
      });
      da.authorise('alice');
    });
  });

  it('Authorisation failed - no auth for one of the files', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var files = [
        {"filepath_by_host" : {"*" : "irods:/seq/foo1"},
          "access_control_group_id": "6"},
        {"filepath_by_host" : {"*" : "irods:/seq/foo2"},
          "access_control_group_id": "8"}                  
                  ];
      var da = new DataAccess(db, files);
      da.on('failed', (username, reason) => {
         expect(username).toBe('alice');
         expect(reason).toBe('Not authorised for some of the files');
         done();
      });
      da.authorise('alice');
    });
  });

  it('Authorisation failed - no auth for any of the files', function(done) {
    MongoClient.connect(url, function(err, db) {
      assert.equal(err, null);
      var files = [
        {"filepath_by_host" : {"*" : "irods:/seq/foo1"},
          "access_control_group_id": "9"},
        {"filepath_by_host" : {"*" : "irods:/seq/foo2"},
          "access_control_group_id": "8"}
                  ];
      var da = new DataAccess(db, files);
      da.on('failed', (username, reason) => {
         expect(username).toBe('alice');
         expect(reason).toBe('Not authorised for any of the files');
         done();
      });
      da.authorise('alice');
    });
  });

  afterAll(() => {
    let out = child.execSync(`mongod --shutdown --dbpath ${tmp_dir} --pidfilepath ${tmp_dir}/mpid`);
    console.log(`\nMONGODB server has been shut down: ${out}`);
    fse.remove(tmp_dir, (err) => {if (err) {console.log(`Error removing ${tmp_dir}: ${err}`)}});
  });
});

