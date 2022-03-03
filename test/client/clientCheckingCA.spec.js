/* globals describe, expect, it, beforeAll, afterAll */

"use strict";

const child  = require('child_process');
const https  = require('https');

const fse = require('fs-extra');
const pem = require('pem');

const config        = require('../../lib/config.js');

let spawn = child.spawn;

describe('test client conneting to server with CA signed cert', () => {
  let ca_cert;
  let server_cert;

  let tmpDir;
  let ca_file_name = 'ca.pem';

  let BASE_PORT  = 5000;
  let PORT_RANGE = 200;
  let SERV_PORT  = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;

  let server;

  beforeAll( (done) => {
    tmpDir = config.tempFilePath('npg_ranger_bin_client_test_');
    fse.ensureDirSync(tmpDir);

    pem.createCertificate({ // For CA
      days:       1,
      commonName: 'CA Certificate'
    }, ( error, cert1) => {
      if (error) {
        done.fail(error);
      }
      ca_cert = cert1;

      fse.writeFileSync( `${tmpDir}/${ca_file_name}`, ca_cert.certificate );

      let certInfo = {
        commonName: 'localhost',
        days:1,
        serviceKey:         ca_cert.serviceKey,
        serviceCertificate: ca_cert.certificate,
        serial: Date.now(),
      };

      pem.createCertificate(certInfo, function (error, cert2) { // For server
        if ( error ) {
          done.fail(error);
        }

        server_cert = cert2;

        server = https.createServer({
          key:  server_cert.clientKey,
          cert: server_cert.certificate,
        }, ( req, res ) => {
          res.end('hi');
        }).listen(SERV_PORT);

        done();
      });
    });
  });

  afterAll( () => {
    server.close();
    try {
      fse.removeSync(tmpDir);
    } catch (e) {
      console.log(e);
    }
  });

  describe('just checking setup', () => {
    it('produces all certs and chain is valid', (done) => {
      expect(ca_cert).toBeDefined();
      expect(server_cert).toBeDefined();

      // Verify the certificate chain is correct
      pem.verifySigningChain( server_cert.certificate, ca_cert.certificate, (error, valid) => {
        if ( error ) {
          done.fail(error);
        }
        expect(valid).toBe(true);
        done();
      });
    });
  });

  describe('test client can connect with CA cert provided', () => {
    it('can connect when providing ca to client', ( done ) => {
      let client = spawn('bin/client.js', [
        `https://localhost:${SERV_PORT}/`,
        `--with_ca=${tmpDir}/${ca_file_name}`
      ]);
      let stdout = '';
      let stderr = '';
      client.stdout.on('data', function(data) {
        stdout += data;
      });
      client.stderr.on('data', function(data) {
        stderr += data;
      });
      client.on('close', function(code) {
        expect(stdout).toEqual('hi');
        expect(stderr).toEqual('');
        expect(code).toBe(0);
        done();
      });
    });
  });

  describe('test client fails connect without CA cert', () => {
    it('fails to connect & shows error message', ( done ) => {
      let client = spawn('bin/client.js', [
        `https://localhost:${SERV_PORT}/`
      ]);
      let stdout = '';
      let stderr = '';
      client.stdout.on('data', function(data) {
        stdout += data;
      });
      client.stderr.on('data', function(data) {
        stderr += data;
      });
      client.on('close', function(code) {
        expect(stdout).toEqual('');
        expect(stderr).toMatch(/unable to verify the first certificate/i);
        expect(code).toBe(1);
        done();
      });
    });
  });
});
