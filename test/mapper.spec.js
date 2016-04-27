"use strict";

const assert      = require('assert');
const child       = require('child_process');
const MongoClient = require('mongodb').MongoClient;
const tmp         = require('tmp');
const fse         = require('fs-extra');

const DataMapper  = require('../lib/mapper.js');

jasmine.DEFAULT_TIMEOUT_INTERVAL = 5000;

const BASE_PORT  = 1100;
const PORT_RANGE = 200;
const PORT = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;
const FIXTURES = 'test/data/fixtures/access_control_group.json';

describe('Data info retrieval', function() {
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

  it('Input validation', function() {
    expect( () => {new DataMapper()} ).toThrowError(ReferenceError,
      'Database handle is required');
    expect( () => {
      let dm = new DataMapper({});
      dm.getFileInfo();
    }).toThrowError(ReferenceError, 'Query object is required');
    expect( () => {
      let dm = new DataMapper({});
      dm.getFileInfo({});
    }).toThrowError(ReferenceError, 'Host name is required');
    expect( () => {
      let dm = new DataMapper({});
      dm.getFileInfo({fruit: "apple"}, 'localhost');
    }).toThrowError(Error,
      'Sample accession number or file name should be given');
  });

  afterAll(() => {
    let out = child.execSync(`mongod --shutdown --dbpath ${tmp_dir} --pidfilepath ${tmp_dir}/mpid`);
    console.log(`\nMONGODB server has been shut down: ${out}`);
    fse.remove(tmp_dir, (err) => {if (err) {console.log(`Error removing ${tmp_dir}: ${err}`)}});
  });
});

