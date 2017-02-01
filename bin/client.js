#!/usr/bin/env node
"use strict";

const assert = require('assert');
const fs     = require('fs');

const cline   = require('commander');
const async   = require("async");
const request = require('request');

const LOGGER        = require('../lib/logsetup.js');
const rangerRequest = require('../lib/client/rangerRequest');
const trailer       = require('../lib/server/http/trailer.js');
const uriUtils      = require('../lib/client/uriUtils.js');

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
  .option('--loglevel <level>', 'level of logging output', /^(error|warn|info|debug)$/i, 'error')
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

var acceptTrailers = cline.acceptTrailers;

LOGGER.level = cline.loglevel;

var url = cline.args[0];
var output;
if ( cline.args.length === 2 ) {
  output = fs.createWriteStream(cline.args[1], {
    flags:     'w',
    autoClose: true
  });
} else {
  output = process.stdout;
}

var exitWithError = (message) => {
  LOGGER.error( message );
  process.exit( 1 );
};

output.on('error', ( err ) => {
  exitWithError( err );
});

output.on('close', () => {
  LOGGER.debug('Output stream has been closed. Exiting...');
  process.exit();
});

const RE_DATA_URI = /^data:/i;

var requestWorker = ( task, callback ) => {
  assert( task.uri, 'uri is required' );
  assert( typeof callback === 'function', 'callback must be of type <function>');

  if ( RE_DATA_URI.test( task.uri ) ) {
    LOGGER.debug('Processing data URI');
    try {
      let buffer = uriUtils.procDataURI( task.uri );
      output.write( buffer );
      callback();
    } catch ( err ) {
      callback( err );
    }
  } else {
    let options = {
      uri:    task.uri,
      method: 'GET'
    };
    options.headers = task.headers ? task.headers : {};
    if ( acceptTrailers ) {
      options.headers.TE = 'trailers';
    }
    let req = request(options);
    req.on('error', ( err ) => {
      LOGGER.error('Error on request ' + err);
      callback( err );
    }).on('response', ( res ) => {
      res.on('error', ( err ) => {
        callback( err );
      });
      if ( res.statusCode === 200 || res.statusCode === 206 ) {
        LOGGER.debug('Status code for <' + task.uri + '>: ' + res.statusCode);
        let contentType = res.headers['content-type'];
        contentType = contentType && contentType.toLowerCase ? contentType.toLowerCase()
                                                             : '';
        if ( contentType.startsWith('application/json') ) {
          try {
            let body = '';
            res.on('data', (data) => {
              body += data;
            });
            res.on('end', () => {
              let uriData = rangerRequest.procJSON( body );
              let q = async.queue( requestWorker, 1 );

              q.drain = () => {
                LOGGER.debug('All items have been processed in internal queue');
              };

              /* jshint -W083 */
              // functions within a loop
              for ( var i = 0; i < uriData.uris.length; i++ ) {
                let newTask = {
                  uri:     uriData.uris[i],
                  headers: uriData.headers4uris[i]
                };
                LOGGER.debug('Pushing to queue: ' + JSON.stringify( newTask ));
                q.push( newTask, ( err ) => {
                  if ( !err ) {
                    LOGGER.debug('Finished task: ' + JSON.stringify( newTask ));
                  } else {
                    callback( err );
                  }
                });
              }
              /* jshint +W083 */
            });
          } catch ( e ) {
            callback( e );
          }
        } else {
          res.on('end', () => {
            LOGGER.debug('End of stream, all data processed');
            if ( acceptTrailers ) {
              LOGGER.debug('Checking trailers');
              let trailerString = trailer.asString(res);
              if (trailerString) {
                LOGGER.info('TRAILERS from ' + task.uri + ': ' + trailerString);
              }
              let dataOK = !trailer.isDataTruncated(options.headers, res);
              if (!dataOK) {
                LOGGER.error('Trailer marked as truncated data. Exiting...');
                callback('Incomplete or truncated data');
              }
            }
            callback();
          });
          res.pipe(output);
        }
      } else {
        callback('Non 200 status ' + res.statusCode);
      }
    });
  }
};

process.nextTick(() => {
  requestWorker({ uri: url }, ( err ) => {
    if ( err ) {
      exitWithError( err );
    } else {
      LOGGER.debug('Success');
      process.exit(0);
    }
  });
});
