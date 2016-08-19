#!/usr/bin/env node
"use strict";

const fs    = require('fs');
const cline = require('commander');

const LOGGER        = require('../lib/logsetup.js');
const RangerRequest = require('../lib/client/rangerRequest');

/**
 * @external fs
 * @see      {@link https://nodejs.org/dist/latest-v4.x/docs/api/fs.html|fs}
 */

/**
 * @external commander
 * @see      {@link https://www.npmjs.com/package/commander|commander}
 */

/**
 * Command line client
 * @module client
 *
 * @requires {@link external:fs|fs}
 * @requires {@link external:commander|commander}
 * @requires module:logsetup
 *
 * Provides a command line client for data retrieval base on GA4GH data sharing
 * API. Is implemented with parallel, asynchronous requests.
 *
 * Data is written to stdout
 *
 * $ client.js "http://192.168.0.1:5050/resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM" > data.bam
 *
 * Or writen to output file if filename is provided as second parameter
 *
 * $ client.js "http://192.168.0.1:5050/resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM" AA0011.bam
 *
 * @author Jaime Tovar
 * @copyright Genome Research Limited 2016
 */

cline
  .version('0.2.2')
  .description('Command line client for GA4GH data streaming')
  .arguments('<url> [output]')
  .option('--debug', 'Show debug output')
  .parse(process.argv);

cline.on('--help', () => {
  console.log('  Examples:');
  console.log('');
  console.log('    $ client.js "http://192.168.0.1:5050/' +
              'resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM"');
  console.log('');
  console.log('    $ client.js "http://192.168.0.1:5050/' +
              'resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM"' +
              ' AA0011.bam');
  console.log('');
});

if ( !cline.args.length ||
     ( cline.args.length != 1 && cline.args.length != 2 ) ) { cline.help(); }

if ( !cline.debug ) {
  LOGGER.level = 'warn';
} else {
  LOGGER.level = 'debug';
}

var url    = cline.args[0];
var output = cline.args.length == 2 ? cline.args[1] : undefined;
var req    = new RangerRequest();

LOGGER.debug('Preparing call');
req.open('GET', url);

req.onreadystatechange = () => {
  LOGGER.debug('readystatechange: ' + req.readyState);
  if ( req.readyState === 4 ) {
    LOGGER.info('Request done with status ' + req.status);
    if ( req.status == 200 || req.status == 206 ) {
      LOGGER.info('Got ' + req.response.byteLength + ' bytes');
      if ( output ) {
        LOGGER.debug('Will write to ' + output);
        fs.open(output, 'w', (err, fd) => {
          if ( err ) {
            LOGGER.error(err);
          }
          fs.write(fd, req.response, 0, req.response.byteLength, (err, written) => { // (err, written, buffer)
            if ( err ) {
              LOGGER.error(err);
            } else {
              LOGGER.debug('Wrote: ' + written + ' bytes to file.');
              fs.close(fd, (err) => {
                LOGGER.error('Error while closing output file ' + err);
              });
            }
          });
        });
      } else {
        LOGGER.debug('Will write to stdout');
        let writeComplete = process.stdout.write(req.response, () => {
          LOGGER.info('Finished writing operation');
        });
        if ( writeComplete ) {
          LOGGER.info('Data was handled completely');
        } else {
          LOGGER.info('Data was not handled completely');
        }
      }
    } else {
      try {
        LOGGER.debug('request: ' + JSON.stringify(req));
      } finally {
        LOGGER.error(req.status + ' ' + req.statusMessage);
      }
    }
  }
};

req.onerror = ( error ) => {
  LOGGER.error('onerror: ' + error);
};

function sendRequest() {
  LOGGER.info('Sending request');
  req.send('');
}

process.nextTick(sendRequest);
