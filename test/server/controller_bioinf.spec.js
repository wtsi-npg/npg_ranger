/* globals describe, it, expect, beforeAll, afterAll */

'use strict';

const assert      = require('assert');
const child       = require('child_process');
const crypto      = require('crypto');
const fse         = require('fs-extra');
const http        = require('http');
const path        = require('path');
const tmp         = require('tmp');
const MongoClient = require('mongodb').MongoClient;

const RangerController = require('../../lib/server/controller.js');
const config           = require('../../lib/config.js');

const BASE_PORT  = 1400;
const PORT_RANGE = 200;
const PORT       = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;
const FIXTURES   = 'test/server/data/fixtures/fileinfo.json';

var tmpDir   = config.tempFilePath('npg_ranger_controller_bioinf_test_');
let db_name  = 'imetacache';
let mongourl = `mongodb://localhost:${PORT}/${db_name}`;
let options  = config.provide( () => { return {mongourl: mongourl}; });

describe('server fetching', () => {
  const server = http.createServer();
  let socket = tmp.tmpNameSync();
  let expectedMd5s = {
    'single file': {
      'SAM' : ['16b3d79daec1da26d98a4e1b63e800b0'],
      'BAM' : ['16b3d79daec1da26d98a4e1b63e800b0'],
      'CRAM': ['16b3d79daec1da26d98a4e1b63e800b0'],
      'VCF' : ['832d2abee762681e6f025ca0df1f38ad']
    },
    'multiple (merged) files': {
      'SAM' : ['79cb05e3fe428da52da346e7d4f6324a', '9b123c8f3a3e8a59584c2193976d1226'],
      'BAM' : ['79cb05e3fe428da52da346e7d4f6324a', '9b123c8f3a3e8a59584c2193976d1226'],
      'CRAM': ['79cb05e3fe428da52da346e7d4f6324a', '9b123c8f3a3e8a59584c2193976d1226'],
      'VCF' : ['832d2abee762681e6f025ca0df1f38ad']
    }
  };
  let testRuns = {
    'single file':             '/file?name=20818_1%23888.bam&format=',
    'multiple (merged) files': '/sample?accession=ABC123456&format='
  };

  beforeAll( (done) => {
    // Start mongod
    fse.ensureDirSync(tmpDir);
    let command = `mongod -f test/server/data/mongodb_conf.yml --port ${PORT} --dbpath ${tmpDir} --pidfilepath ${tmpDir}/mpid --logpath ${tmpDir}/dbserver.log`;
    console.log(`\nCommand to start MONGO DB daemon: ${command}`);
    let out = child.execSync(command);
    console.log(`Started MONGO DB daemon: ${out}`);

    // Import data set from fixtures.json to fileinfo
    command = `mongoimport --port ${PORT} --db ${db_name} --collection fileinfo --jsonArray --file ${FIXTURES}`;
    out = child.execSync(command);
    console.log(`Loaded data to MONGO DB: ${out}`);

    // Uncomment these lines to view server debug output
    //let LOGGER = require('../../lib/logsetup.js');
    //LOGGER.level = 'debug';

    MongoClient.connect(mongourl, (err, db) => {
      assert.equal(err, null, `failed to connect to ${mongourl}: ${err}`);

      // The model runs samtools merge in a temporary directory,
      // so mongo needs to return an absolute path to the data.
      // Difficult because tests need to be portable.
      // So, test script updates the entries with current directory.
      let cwd = process.cwd();
      let collection = db.collection('fileinfo');

      let updatePromises = ['20818_1#888.bam', '20907_1#888.bam']
        .map(function createUpdatePromise(dataObj) {
          return collection.findOne({'data_object': dataObj})
          .then(function docFound(doc) {
            collection.findOneAndUpdate(
              {'data_object': dataObj},
              {'$set':
                { 'filepath_by_host.*': path.join(cwd, doc.filepath_by_host['*']) }
              }
            );
          }, function docNotFound(reason) {
            console.log('Document was not found: ' + reason);
          });
        });

      let dbClose = (dbConn) => {
        if (dbConn) {
          try {
            dbConn.close();
          } catch (err) { console.log(err); }
        }
      };

      server.on('close', () => {
        dbClose(db);
      });

      let listenPromise = new Promise( (resolve) => {
        server.listen(socket, () => {
          console.log('server listening on ' + socket);
          resolve();
        });
      });

      server.on('request', (request, response) => {
        let c = new RangerController(request, response, db, tmpDir, true);
        c.handleRequest(options.get('hostname'));
      });

      Promise.all(updatePromises.concat([listenPromise]))
      .then(function serverReady() {
        done();
      }, ( reason ) => {
        console.log('Server wasn\'t ready: ' + reason);
      });
    });
  });

  afterAll( () => {
    server.close( function closed(err){
      if (err) {
        console.log('Tried to close server, but it was already closed.');
      }
      child.execSync(`mongo 'mongodb://localhost:${PORT}/admin' --eval 'db.shutdownServer()'`);
      // the above shutdown command can return before server is
      // shut down, so set a short timeout to make sure
      setTimeout( () => {
        try { fse.unlinkSync(socket); } catch (e) { console.log(e); }
        fse.removeSync(tmpDir);
      }, 1000);
    });
  });

  Object.keys(testRuns).forEach( ( description ) => {
    describe(description, () => {
      Object.keys(expectedMd5s[description]).forEach( ( format ) => {
        it('run controller on ' + description + ' outputting ' + format, ( done ) => {
          function isOneOf(subject, expecteds) {
            let matched = false;
            for (let i = 0; i < expecteds.length; i++) {
              if (subject === expecteds[i]) {
                matched = true;
                break;
              }
            }
            return matched;
          }
          http.get(
              {
                socketPath: socket,
                path: testRuns[description] + format,
                headers: {TE: 'trailers'}
              }, ( res ) => {

            let hash = crypto.createHash('md5');

            if (format !== 'VCF') {
              // The header will change between runs, so use bamseqchksum.
              // Bamseqchksum is hard to compare, so md5 it and compare that.
              let bamseqchksum = child.spawn('bamseqchksum',
                [
                  'inputformat=' + format.toLowerCase(),
                  'reference=test/server/data/references/PhiX/all/fasta/phix_unsnipped_short_no_N.fa'
                ]);
              res.on('data', ( data ) => {
                bamseqchksum.stdin.write(data);
              });
              bamseqchksum.stdout.on('data', ( data ) => {
                hash.update(data);
              });
              res.on('end', () => {
                bamseqchksum.stdin.end();
              });
              bamseqchksum.stdout.on('end', () => {
                let hashDigest = hash.digest('hex');
                let match = isOneOf(hashDigest, expectedMd5s[description][format]);
                expect(match).toBe(true);
                done();
              });
            } else {
              let body = '';
              let file = fse.createWriteStream('./out.vcf');
              res.on('data', ( data ) => {
                body += data;
              });
              res.on('end', () => {
                // in VCF file, header is unpredictable, so remove and md5 remaining data
                file.end(body);
                let hashDigest = hash.update(body.replace(/#.*?\n/g, ''))
                                     .digest('hex');
                let match = isOneOf(hashDigest,
                                    expectedMd5s[description][format]
                                   );
                expect(match).toBe(true);
                done();
              });
            }
          });
        });
      });
    });
  });
});
