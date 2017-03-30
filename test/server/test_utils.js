"use strict";

const fs  = require('fs');
const fse = require('fs-extra');
const pem = require('pem');

const KEY_EXT  = 'key';
const CERT_EXT = 'cert';

/**
 * @module test/server/test_utils
 *
 * @description
 * <p>Common functions for test modules</p>
 *
 * @requires {@link external:fs|fs}
 * @requires {@link external:fs-extra|fs-extra}
 * @requires {@link external:pem|pem}
 *
 * @author Jaime Tovar
 * @copyright Genome Research Limited 2017
 */

/**
 * Best-effort approach to remove socket files from the system. If something
 * goes wrong it will report the error to console.
 * @param  {string} socket socket path
 */
var removeSocket = ( socket ) => {
  try {
    // check if constants is defined for compatibility between node v4 and v6.3+
    let flag = fs.constants ? fs.constants.W_OK : fs.W_OK;
    fs.access(socket, flag, ( err ) => {
      if ( !err ) {
        fs.unlinkSync(socket);
      } else {
        console.log( err );
      }
    });
  } catch (e) { console.log(e); }
};

/**
 * Generates a certificate chain with one CA and two certificates signed by that
 * CA. Both a certificate and a key will be generated for the children
 * certificates. Certificates will have the ".cert" filename extension, keys
 * will have the ".key" filename extension.
 * @param  {String}   path         where files will be saved
 * @param  {String}   ca_prefix    prefix for the CA filename
 * @param  {String}   cert1_prefix prefix for the 1st certificate filenames
 * @param  {String}   cert2_prefix prefix for the 2nd certficate filenames
 * @param  {Function} callback     will be called with the cerficate objects
 *                                 with (err, CA, cert1, cert2) signature.
 */
let create_certificates = (path, ca_prefix, cert1_prefix, cert2_prefix, callback) => {
  fse.ensureDirSync(path);

  let ca_file_name = `${ca_prefix}.${CERT_EXT}`;

  pem.createCertificate({ // For CA
    days:       1,
    commonName: 'CA Certificate'
  }, ( error, cert0) => {
    if (error) {
      callback(error);
    }
    let ca_cert = cert0;

    try {
      fse.writeFileSync( `${path}/${ca_file_name}`, ca_cert.certificate );
    } catch ( e ) {
      callback( e );
    }

    let certInfo = {
      commonName: 'localhost',
      days:1,
      serviceKey:         ca_cert.serviceKey,
      serviceCertificate: ca_cert.certificate,
      serial: Date.now(),
    };

    pem.createCertificate(certInfo, function (error, cert1) { // For cert 1
      if ( error ) {
        callback(error);
      }

      try {
        fse.writeFileSync(`${path}/${cert1_prefix}.${CERT_EXT}`, cert1.certificate);
        fse.writeFileSync(`${path}/${cert1_prefix}.${KEY_EXT}`, cert1.clientKey);
      } catch ( e ) {
        callback( e );
      }

      let certInfo2 = {
        commonName: 'localhost',
        days:1,
        serviceKey:         ca_cert.serviceKey,
        serviceCertificate: ca_cert.certificate,
        serial: Date.now(),
      };

      pem.createCertificate(certInfo2, function (error, cert2) { // For cert 2
        if ( error ) {
          callback(error);
        }

        try {
          fse.writeFileSync(`${path}/${cert2_prefix}.${CERT_EXT}`, cert2.certificate);
          fse.writeFileSync(`${path}/${cert2_prefix}.${KEY_EXT}`, cert2.clientKey);
        } catch ( e ) {
          callback( e );
        }

        callback(null, cert0, cert1, cert2);
      });
    });
  });
};

module.exports = {
  create_certificates: create_certificates,
  removeSocket: removeSocket
};
