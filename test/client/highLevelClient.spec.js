/* globals describe, xdescribe, expect, it, fail, beforeAll, afterAll */

"use strict";

const assert = require('assert');
const child  = require('child_process');
const crypto = require('crypto');
const fse    = require('fs-extra');
const md5    = require('js-md5');
const path   = require('path');
const MongoClient = require('mongodb').MongoClient;

const config        = require('../../lib/config.js');
const RangerRequest = require('../../lib/client/rangerRequest');

// TODO remove these tests when finished building new
xdescribe('Testing high level', () => {
  it('Success with Google', ( done ) => {
    var req = new RangerRequest();
    // Google
    var url = 'http://104.196.18.135/readgroupsets/CMvnhpKTFhD04eLE-q2yxnU?referenceName=1&start=167856&end=173507&format=BAM';

    var expected = {
      bytes: 63778,
      md5: '94a3fc74146898aab618cad666a3a54a'
    };

    req.open('GET', url);

    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.status === 200 || req.status === 206).toBe(true);
        expect(req.response.byteLength).toBe(expected.bytes, 'Got correct number of bytes');
        expect(md5(new Uint8Array(req.response))).toBe(expected.md5, 'Data hash matches expected');
        done();
      }
    };

    req.send('');
  }, 5000);

  it('Success with Google, but with different data', ( done ) => {
    var req = new RangerRequest();
    // Google
    var url = 'http://104.196.18.135/readgroupsets/CMvnhpKTFhD04eLE-q2yxnU?referenceName=1&start=160000&end=165000&format=BAM';

    var expected = {
      bytes: 57977,
      md5: '41a7756429c37a78b65b5aa4a5891e54'
    };

    req.open('GET', url);

    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.status === 200 || req.status === 206).toBe(true);
        expect(req.response.byteLength).toBe(expected.bytes, 'Got correct number of bytes');
        expect(md5(new Uint8Array(req.response))).toBe(expected.md5, 'Data hash matches expected');
        done();
      }
    };

    req.send('');
  }, 5000);

  it('Success with DNANexus - 1k genomes', ( done ) => {
    var req = new RangerRequest();
    // DNANexus
    var url = 'http://htsnexus.rnd.dnanex.us/v1/reads/1000genomes_low_coverage/NA20276?referenceName=22&start=16100000&end=16105000&format=BAM';

    var expected = {
      bytes: 279561,
      md5: 'd7924c94113ca4e76d64d321f04e6766'
    };

    req.open('GET', url);

    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.status === 200 || req.status === 206).toBe(true);
        expect(req.response.byteLength).toBe(expected.bytes, 'Got correct number of bytes');
        expect(md5(new Uint8Array(req.response))).toBe(expected.md5, 'Data hash matches expected');
        done();
      }
    };

    req.send('');
  }, 5000);

  it('Success with DNANexus - platinum', ( done ) => {
    var req = new RangerRequest();
    // DNANexus
    var url = 'http://htsnexus.rnd.dnanex.us/v1/reads/platinum/NA12891?referenceName=chr22&start=16100000&end=16105000&format=BAM';

    var expected = {
      bytes: 814155,
      md5: 'ca298fbb2346bd2df88e647d2e739827'
    };

    req.open('GET', url);

    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.status === 200 || req.status === 206).toBe(true);
        expect(req.response.byteLength).toBe(expected.bytes, 'Got correct number of bytes');
        expect(md5(new Uint8Array(req.response))).toBe(expected.md5, 'Data hash matches expected');
        done();
      }
    };

    req.send('');
  }, 20000);

  it('Not success with wrong path', ( done ) => {
    var req = new RangerRequest();
    // Google
    var url = 'http://104.196.18.135/readgroupset/CMvnhpKTFhD04eLE-q2yxnU?referenceName=1&start=167856&end=173507&format=BAM';

    req.open('GET', url);

    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.status).toBe(404, 'Correct status code');
        done();
      }
    };

    req.send('');
  }, 5000);

  it('Not success with wrong server', ( done ) => {
    var req = new RangerRequest();
    // Google
    var url = 'http://127.0.0./someData';

    req.open('GET', url);

    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.status).toBe(424, 'Correct status code for server address error');
        expect(req.statusMessage).toMatch(/ENOTFOUND/, 'Error message includes statusMessage');
        expect(req.statusMessage).toMatch(/127\.0\.0\./, 'Error message includes the url requested');
        done();
      }
    };

    req.send('');
  }, 5000);
});

describe('Running with ranger server', () => {
  let spawn    = child.spawn;
  let execSync = child.execSync;

  let BASE_PORT  = 1400;
  let PORT_RANGE = 200;
  let SERV_PORT  = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;
  let MONGO_PORT = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;
  const FIXTURES = 'test/server/data/fixtures/fileinfo.json';

  let dbName        = 'imetacache';
  let mongourl      = `mongodb://localhost:${MONGO_PORT}/${dbName}`;
  let tmpDir        = config.tempFilePath('npg_ranger_bin_server_test_');
  let serverCommand = 'bin/server.js';

  beforeAll((done) => {
    // Start mongod
    fse.ensureDirSync(tmpDir);
    let command = `mongod -f test/server/data/mongodb_conf.yml --port ${MONGO_PORT} --dbpath ${tmpDir} --pidfilepath ${tmpDir}/mpid --logpath ${tmpDir}/dbserver.log`;
    console.log(`\nCommand to start MONGODB daemon: ${command}`);
    let out = execSync(command);
    console.log(`Started MONGODB daemon: ${out}`);

    // Import data set from fixtures.json
    command = `mongoimport --port ${MONGO_PORT} --db ${dbName} --collection fileinfo --jsonArray --file ${FIXTURES}`;
    out = execSync(command);
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
        fail('Didn\'t update all docs: ' + reason);
        done();
      });
    });
  });

  afterAll( (done) => {
    child.execSync(`mongo '${mongourl}' --eval 'db.shutdownServer()'`);
    setTimeout( () => {
      fse.removeSync(tmpDir);
      done();
    }, 1000);
  });

  it('Running with ranger client', (done) => {
    let serv = spawn(serverCommand, [
      '-s',
      '-d',
      '-n0',
      `-p${SERV_PORT}`,
      `${mongourl}`]);
    serv.on('close', (code, signal) => {
      if (code || signal) {
        console.log(code || signal);
      }
      done();
    });

    serv.stdout.on('data', (data) => {
      if (data.toString().match(/Server listening on /)) {
        // Server is listening and ready for connection
        let client = spawn('bin/client.js', [
          `http://localhost:${SERV_PORT}/file?name=20818_1%23888.bam&format=SAM`]);
        let bamseqchksum = spawn('bamseqchksum', ['inputformat=sam']);
        client.stdout.pipe(bamseqchksum.stdin);
        let hash = crypto.createHash('md5');
        bamseqchksum.stdout.on('data', (data) => {
          hash.update(data.toString());
        });
        bamseqchksum.on('exit', () => {
          expect(hash.digest('hex')).toBe('16b3d79daec1da26d98a4e1b63e800b0');
          serv.kill();
        });
      }
    });
  }, 20000);
});
