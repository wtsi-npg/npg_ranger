/* globals describe, it, expect, beforeAll, afterAll */

'use strict';

const assert      = require('assert');
const child       = require('child_process');
const crypto      = require('crypto');
const fse         = require('fs-extra');
const http        = require('http');
const MongoClient = require('mongodb').MongoClient;
const path        = require('path');
const url         = require('url');

const config      = require('../../lib/config.js');

const BASE_PORT  = 1400;
const PORT_RANGE = 200;
const MONGO_PORT = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;
const SERV_PORT  = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;
const FIXTURES   = 'test/server/data/fixtures/fileinfo.json';

let tmpDir   = config.tempFilePath('npg_ranger_bin_server_test_');
let db_name  = 'imetacache';
let mongourl = `mongodb://localhost:${MONGO_PORT}/${db_name}`;

describe('Running server bin', () => {
  beforeAll( (done) => {
    // Start mongod
    fse.ensureDirSync(tmpDir);
    let command = `mongod -f test/server/data/mongodb_conf.yml --port ${MONGO_PORT} --dbpath ${tmpDir} --pidfilepath ${tmpDir}/mpid --logpath ${tmpDir}/dbserver.log`;
    console.log(`\nCommand to start MONGODB daemon: ${command}`);
    let out = child.execSync(command);
    console.log(`Started MONGODB daemon: ${out}`);

    // Import data set from fixtures.json
    command = `mongoimport --port ${MONGO_PORT} --db ${db_name} --collection fileinfo --jsonArray --file ${FIXTURES}`;
    out = child.execSync(command);
    console.log(`Loaded data to MONGO DB: ${out}`);

    MongoClient.connect(mongourl, (err, db) => {
      assert.equal(err, null, `failed to connect to ${mongourl}: ${err}`);

      // The model runs samtools merge in a temporary directory,
      // so mongo needs to return an absolute path to the data.
      // Difficult because tests need to be portable.
      // So, test script updates the entries with current directory.
      let cwd = process.cwd();
      let collection = db.collection('fileinfo');

      let updatePromises = ['20818_1#888.bam', '20907_1#888.bam']
        .map(function(dataObj) {
          return collection.findOne({'data_object': dataObj})
          .then(function(doc) {
            collection.findOneAndUpdate(
              {'data_object': dataObj},
              {'$set':
                {'filepath_by_host.*': path.join(cwd, doc.filepath_by_host['*'])}
              }
            );
          }, function MongoDocNotFound(reason) {
            console.log('Document for ' + dataObj + ' was not found: ' + reason);
          });
        });

      Promise.all(updatePromises).then( () => {
        done();
      }, (reason) => {
        console.log('Didn\'t update all docs: ' + reason);
      });
    });
  });

  afterAll( () => {
    child.execSync(`mongo 'mongodb://localhost:${MONGO_PORT}/admin' --eval 'db.shutdownServer()'`);
    setTimeout( () => {
      fse.removeSync(tmpDir);
    }, 1000);
  });


  it('run the server', (done) => {
    let serv = child.fork('bin/server.js', [`-sp${SERV_PORT}`, '-m',`mongodb://localhost:${MONGO_PORT}/imetacache`]);
    serv.on('close', (code, signal) => {
      if (code || signal) {
        console.log(code || signal);
      }
      done();
    });
    serv.on('message', (m) => {
      if (m.listening === true) {
        http.get({path: '/file?name=20818_1%23888.bam&format=SAM', port: SERV_PORT}, (res) => {
          let bamseqchksum = child.spawn('bamseqchksum', ['inputformat=sam']);
          res.pipe(bamseqchksum.stdin);
          let hash = crypto.createHash('md5');
          bamseqchksum.stdout.on('data', (data) => {
            hash.update(data.toString());
          });
          bamseqchksum.on('exit', () => {
            expect(hash.digest('hex')).toBe('16b3d79daec1da26d98a4e1b63e800b0');
            serv.kill();
          });
        });
      }
    });
  });

  it('run the server, connect via socket', (done) => {
    let serv = child.fork('bin/server.js', ['-sp', `${tmpDir}/npg_ranger.sock`, '-m', `mongodb://localhost:${MONGO_PORT}/imetacache`]);
    serv.on('close', (code, signal) => {
      if (code || signal) {
        console.log(code || signal);
      }
      done();
    });
    serv.on('message', (m) => {
      if (m.listening === true) {
        http.get({path: '/file?name=20818_1%23888.bam&format=SAM', socketPath: `${tmpDir}/npg_ranger.sock`}, (res) => {
          let bamseqchksum = child.spawn('bamseqchksum', ['inputformat=sam']);
          res.pipe(bamseqchksum.stdin);
          let hash = crypto.createHash('md5');
          bamseqchksum.stdout.on('data', (data) => {
            hash.update(data.toString());
          });
          bamseqchksum.on('exit', () => {
            expect(hash.digest('hex')).toBe('16b3d79daec1da26d98a4e1b63e800b0');
            serv.kill();
          });
        });
      }
    });
  });

  it('Run server, merge files', (done) => {
    let serv = child.fork('bin/server.js', [`-sp${SERV_PORT}`, '-m', `mongodb://localhost:${MONGO_PORT}/imetacache`]);
    serv.on('close', (code, signal) => {
      if (code || signal) {
        console.log(code || signal);
      }
      done();
    });
    serv.on('message', (m) => {
      if (m.listening === true) {
        http.get({path: '/sample?accession=ABC123456&format=SAM', port: SERV_PORT}, (res) => {
          let bamseqchksum = child.spawn('bamseqchksum', ['inputformat=sam']);
          res.pipe(bamseqchksum.stdin);
          let hash = crypto.createHash('md5');
          bamseqchksum.stdout.on('data', (data) => {
            hash.update(data.toString());
          });
          bamseqchksum.on('exit', () => {
            let chksums = ['79cb05e3fe428da52da346e7d4f6324a', '9b123c8f3a3e8a59584c2193976d1226'];
            expect(chksums).toContain(hash.digest('hex'));
            serv.kill();
          });
        });
      }
    });
  });

  it('Run server, JSON redirect response', (done) => {
    let serv = child.fork('bin/server.js', [`-sp${SERV_PORT}`, '-m', `mongodb://localhost:${MONGO_PORT}/imetacache`]);
    serv.on('close', (code, signal) => {
      if (code || signal) {
        console.log(code || signal);
      }
      done();
    });
    serv.on('message', (m) => {
      if (m.listening === true) {
        http.get({path: '/ga4gh/v.0.1/get/sample/ABC123456', port: SERV_PORT}, (res) => {
          let jsonResponse = '';
          res.on('data', (data) => {
            jsonResponse += data.toString();
          });
          res.on('end', () => {
            let responseObj;
            expect( () => {responseObj = JSON.parse(jsonResponse);}).not.toThrow();
            expect(responseObj.format).toBe('BAM');
            expect(url.parse(responseObj.urls[0].url).path).toBe('/sample?accession=ABC123456&format=BAM');
            serv.kill();
          });
        });
      }
    });
  });
});
