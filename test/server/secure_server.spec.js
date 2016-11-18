/* globals describe, it, expect, beforeAll, afterAll, afterEach, fail */

"use strict";

const child  = require('child_process');
const https  = require('https');

const fse = require('fs-extra');
const pem = require('pem');

const config = require('../../lib/config.js');

var tmpDir = config.tempFilePath('npg_ranger_ssl_test_');

describe('test running https server', () => {

  const BASE_PORT  = 1400;
  const PORT_RANGE = 100;
  const MONGO_PORT = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;
  const SERV_PORT  = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT + PORT_RANGE;
  const FIXTURES   = 'test/server/data/fixtures/fileinfo_secure_server.json';

  let serverCommand = 'bin/server.js';
  let dbName        = 'npg_ranger_test';
  let mongourl      = `mongodb://localhost:${MONGO_PORT}/${dbName}`;

  let private_pem = tmpDir + '/server_key.pem';
  let cert_pem    = tmpDir + '/server_cert.pem';

  let spawn = child.spawn;

  let serv;

  let startServer = ( myDone ) => {
    let s = spawn(serverCommand, [
      '-s',
      '-d',
      '--startssl',
      '--secure_key', private_pem,
      '--secure_cert', cert_pem,
      '-n0',
      `-p${SERV_PORT}`,
      `-m${mongourl}`]);
    s.on('close', (code, signal) => {
      console.log('Closing server');
      if (code || signal) {
        myDone.fail('Server failed with error: ' + (code || signal));
      }
    });
    return s;
  };

  afterEach( () => {
    serv.kill();
  });

  afterAll( ( done ) => {
    child.execSync(`mongo 'mongodb://localhost:${MONGO_PORT}/admin' --eval 'db.shutdownServer()'`);
    console.log('\nMONGODB server has been shut down');

    setTimeout( () => {
      fse.removeSync(tmpDir);
      done();
    }, 500);
  });

  beforeAll( ( done ) => {
    fse.ensureDirSync(tmpDir);
    console.log(tmpDir);

    let certPromise = new Promise( ( resolve, reject ) => {
      pem.createCertificate( { days:2, selfSigned:true }, ( err, keys ) => {
        if ( err ) {
          reject( err );
        }
        try {
          fse.writeFileSync(private_pem, keys.serviceKey);
        } catch (e) {
          reject(e);
        }
        try {
          fse.writeFileSync(cert_pem, keys.certificate);
        } catch (e) {
          reject(e);
        }
        resolve( keys );
      });
    });

    certPromise.then( () => {
      console.log(`MONGO data directory: ${tmpDir}`);
      let command = `mongod -f test/server/data/mongodb_conf.yml --port ${MONGO_PORT} --dbpath ${tmpDir} --pidfilepath ${tmpDir}/mpid --logpath ${tmpDir}/dbserver.log`;
      console.log(`\nCommand to start MONGO DB daemon: ${command}`);
      let out = child.execSync(command);
      console.log(`Started MONGO DB daemon: ${out}`);
      command = `mongoimport --port ${MONGO_PORT} --db ${dbName} --collection fileinfo --jsonArray --file ${FIXTURES}`;
      out = child.execSync(command);
      console.log(`Loaded data to MONGO DB: ${out}`);

      done();
    }).catch( ( reason ) => {
      fail(reason);
    });
  });

  it('can reply with a reference from a https requests', ( done ) => {
    serv = startServer( done );
    let acc  = 'XYZ120923';
    serv.stdout.on('data', (data) => {
      if (data.toString().match(/Server listening on /)) {
        // Server is listening and ready for connection
        let options = {
          hostname:           'localhost',
          port:               SERV_PORT,
          path:               `/sample/${acc}/reference`,
          method:             'GET',
          rejectUnauthorized: false // so the client disregards self signed cert
        };
        options.agent = new https.Agent(options);

        let req = https.request(options, (res) => {
          expect(res.statusCode).toBe(200);
          try {
            // Check expected errors from a self signed cert
            expect(res.client.authorized).toBe(false);
            expect(res.client.authorizationError).toBe(
              'DEPTH_ZERO_SELF_SIGNED_CERT'
            );
            // Check it was encrypted
            expect(res.client.encrypted).toBe(true);
            expect(res.headers['content-type']).toEqual('application/json');
          } catch (e) {
            fail(e);
          }

          let body ='';
          res.on('data', (d) => {
            try { body += d; } catch (e) { fail(e); }
          });

          res.on('end', () => {
            try {
              let jsonData = JSON.parse(body);
              expect(jsonData.accession).toBeDefined();
              expect(jsonData.reference).toBeDefined();
              expect(jsonData.accession).toEqual(acc);
              expect(jsonData.reference).toEqual(
                '/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'
              );
              done();
            } catch (e) {
              fail(e);
            }
          });
        });
        req.on('error', (e) => {
          fail(e);
        });
        req.end();
      }
    });
  });
});
