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
 * @module client
 *
 * @requires {@link external:fs|fs}
 * @requires {@link external:commander|commander}
 * @requires {@link module:logsetup|logsetup}
 *
 * @description
 * <p>Command line client</p>
 *
 * <p>A command line client for sequencing data retrieval base on GA4GH data sharing
 * API. The client is implemented with parallel, asynchronous requests to every URL
 * listed in the initial JSON response of the GA4GH-compliant server. The clent
 * itself is synchronous, ie the data is returned to the caller only when all requests
 * successfully completed. Untill that moment all received chunks of data are
 * kept in memory, thus limiting the total amount of data this client can retrieve.</p>
 *
 * <p>The client is also able to process direct HTTP GET requests.</p>
 *
 * <p>Data is written to stdout</p>
 *
 * <code>$ client.js "http://some_server_url/resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM" > data.bam</code>
 *
 * <p>... or to an output file if a file path is provided as the second parameter.</p>
 *
 * <code>$ client.js "http://some_server_url/resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM" AA0011.bam</code>
 *
 * <p>If the response was empty, an empty file is created.</p>
 *
 * <p>The client exits with an error code <em>1</em> if an error occured when
 * requesting/receiving the data. If <em>--accept-trailers</em> option is
 * enabled:</p>
 * <ol>
 *   <li>for each response, trailers, if available, are written to
 *       <em>stderr</em></li>
 *   <li>if any of the responses has <em>data-truncated</em> trailer set to
 *       true, the script exits with an error code <em>1</em></li>
 * </ol>
 *
 * @author Marina Gourtovaia
 * @author Jaime Tovar
 * @copyright Genome Research Limited 2016
 */

cline
  .version(require('../package.json').version)
  .description('Command line client for GA4GH data streaming')
  .arguments('<url> [output]')
  .option('--loglevel <level>', 'level of logging output', /^(error|warn|info|debug)$/i, 'warn')
  .option('--accept-trailers', 'Request trailers from server')
  .parse(process.argv);

cline.on('--help', () => {
  console.log('  Examples:');
  console.log('');
  console.log('    $ client.js "http://some_server_url/' +
              'resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM"');
  console.log('');
  console.log('    $ client.js "http://some_server_url/' +
              'resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM"' +
              ' AA0011.bam');
  console.log('');
  console.log('  If you know the server supports trailers, we suggest you execute' +
              ' with "--accept-trailers" option to improve error control.');
  console.log('');
  console.log('    $ client.js --accept-trailers "http://some_server_url/' +
              'resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM"' +
              ' AA0011.bam');
  console.log('');
});

if ( !cline.args.length ||
     ( cline.args.length != 1 && cline.args.length != 2 ) ) { cline.help(); }

LOGGER.setup(true);

LOGGER.level = cline.loglevel;

var url    = cline.args[0];
var output = cline.args.length == 2 ? cline.args[1] : undefined;

LOGGER.debug(`Preparing request for ${url}`);
var req    = new RangerRequest( cline.acceptTrailers );
req.open('GET', url);

var exitWithError = (message) => {
  LOGGER.error(message);
  process.exit(1);
};

req.onreadystatechange = () => {
  LOGGER.debug('readystatechange: ' + req.readyState);
  if ( req.readyState === 4 ) {
    LOGGER.info('Request done with status ' + req.status);
    if ( req.status == 200 || req.status == 206 ) {
      let numBytes = req.response ? req.response.byteLength : 0;
      LOGGER.info(`Got ${numBytes} bytes`);
      if ( output ) {
        LOGGER.debug('Will write to ' + output);
        fs.open(output, 'w', (err, fd) => {
          if ( err ) {
            exitWithError(`Failed to open ${output} for writing: ` + err);
          }
          let content = req.response || '';
          fs.write(fd, content, 0, numBytes, (err, written) => {
            if ( err ) {
              exitWithError(`Failed to write to ${output}: ` + err);
            } else {
              LOGGER.debug('Wrote: ' + written + ' bytes to file.');
              fs.close(fd, (err) => {
                if ( err ) {
                  LOGGER.error('Error while closing output file ' + err);
                } else {
                  LOGGER.debug('Closed output file');
                }
              });
            }
          });
        });
      } else {
        if (numBytes) {
          LOGGER.debug('Will write to stdout');
          process.stdout.on('error', (err) => {
            exitWithError('Failed to write complete data to stdout: ' + err);
          });
          process.stdout.write(req.response, () => {
            LOGGER.info('Finished writing operation');
          });
        }
      }
    } else {
      try {
        LOGGER.debug('request: ' + JSON.stringify(req));
      } finally {
        let s = req.status;
        let m = req.statusMessage;
        exitWithError(`Request failed, status ${s}, message '${m}'`);
      }
    }
  }
};

req.onerror = ( error ) => {
  exitWithError('Error dealing with request: ' + error);
};

process.nextTick(() => { // All callbacks are set now
  LOGGER.info('Sending request');
  req.send('');
});
