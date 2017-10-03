#!/usr/bin/env node
"use strict";

const assert = require('assert');
const fs     = require('fs');
const path   = require('path');

const PassThrough = require('stream').PassThrough;

const cline   = require('commander');
const async   = require("async");
const request = require('request');

const LOGGER        = require('../lib/logsetup.js');
const rangerRequest = require('../lib/client/rangerRequest');
const trailer       = require('../lib/server/http/trailer.js');
const tokenUtils    = require('../lib/token_utils');
const uriUtils      = require('../lib/client/uriUtils.js');
const constants     = require('../lib/constants');

/**
 * @external assert
 * @see      {@link https://nodejs.org/dist/latest-v4.x/docs/api/assert.html|assert}
 */

/**
 * @external fs
 * @see      {@link https://nodejs.org/dist/latest-v4.x/docs/api/fs.html|fs}
 */

/**
 * @external path
 * @see      {@link https://nodejs.org/dist/latest-v4.x/docs/api/path.html|path}
 */

/**
 * @external stream
 * @see      {@link https://nodejs.org/dist/latest-v4.x/docs/api/stream.html|}
 */

/**
 * @external async
 * @see      {@link https://www.npmjs.com/package/async|async}
 */

/**
 * @external request
 * @see      {@link https://www.npmjs.com/package/request|request}
 */

/**
 * @external commander
 * @see      {@link https://www.npmjs.com/package/commander|commander}
 */

/**
 * @module client
 *
 * @requires {@link external:fs|fs}
 * @requires {@link external:path|path}
 * @requires {@link external:stream|stream}
 * @requires {@link external:commander|commander}
 * @requires {@link module:logsetup|logsetup}
 *
 * @description
 * <p>Command line client</p>
 *
 * <p>A command line client for sequencing data retrieval base on GA4GH data sharing
 * API. The client sequentially processes every URL listed in the initial
 * JSON response of the GA4GH-compliant server. The responses are processed and
 * data is streamed to the specified output (file or stdout).</p>
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
 * <p>If a token is required to authorise access to resources, a path to the
 * configuration file in <em>JSON</em> format can be provided with the
 * <em>--token_config</em> option.</p>
 *
 * <code>$ client.js --token_config path/to/file.json "http://some_server..."</code>
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
 * <p>Please take into consideration the current implementation is recursive by
 * design. Therefore, circular references in <em>JSON</em> resources will
 * produce infinite (or very large number of) requests to the resources.</p>
 *
 * @copyright Genome Research Limited 2017
 */

const TOKEN_BEARER_KEY_NAME = constants.TOKEN_BEARER_KEY_NAME;
const TOKEN_CONFIG_KEY_NAME = constants.TOKEN_CONFIG_KEY_NAME;

cline
  .version(require('../package.json').version)
  .description('Command line client for GA4GH data streaming')
  .arguments('<url> [output]')
  .option('--accept-trailers', 'Request trailers from server')
  .option('--loglevel <level>', 'level of logging output', /^(error|warn|info|debug)$/i, 'error')
  .option('--token_config <token_config_file>', 'path to file with token configuration in json format')
  .option('--with_ca <path_to_ca_file>', 'path to CA file')
  .parse(process.argv);

cline.on('--help', () => {
  console.log('  Examples:');
  console.log('');
  console.log('    $ bin/client.js "http://some_server_url/' +
              'resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM"');
  console.log('');
  console.log('    $ bin/client.js "http://some_server_url/' +
              'resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM"' +
              ' AA0011.bam');
  console.log('');
  console.log('  If the server requires a token for authorisation you can use' +
              ' --token_config <file>. To provide a path for the json file storing' +
              ' the configuration.');
  console.log('');
  console.log('    $ bin/client.js --token_config path/to/file.json "http://some_server_url/' +
              'resources/AA0011?referenceName=1&start=500&end=1000&format=BAM"');
  console.log('');
  console.log('  Whenever a token is being used HTTPS should also be used. If the' +
              ' server is configured with a private CA you will need to provide' +
              ' the compatible CA certificate to establish the connection.');
  console.log('');
  console.log('    $ bin/client.js --with_ca path/to/your_ca.crt --token_config' +
              ' path/to/file.json "http://some ..."');
  console.log('');
  console.log('  If you know the server supports trailers, we suggest you execute' +
              ' with "--accept-trailers" option to improve error control.');
  console.log('');
  console.log('    $ bin/client.js --accept-trailers "http://some_server_url/' +
              'resources/AA0011?referenceName=1&start=167856&end=173507&format=BAM"' +
              ' AA0011.bam');
  console.log('');
});

if ( !cline.args.length ||
     ( cline.args.length != 1 && cline.args.length != 2 ) ) { cline.help(); }

var acceptTrailers = cline.acceptTrailers;

var token_config = cline.token_config;
var token;

LOGGER.level = cline.loglevel;

var ca_file = cline.with_ca;

var url = cline.args[0];
var output = new PassThrough();
if ( cline.args.length === 2 ) {
  var fileoutput = fs.createWriteStream(cline.args[1], {
    flags:     'w',
    autoClose: true
  });
  output.pipe(fileoutput);
} else {
  process.stdout.on('error', err => {
    if (err.code == "EPIPE") {
      // next process in the pipe closed e.g. samtools printing only headers
      process.exit(0);
    }
  });
  output.pipe(process.stdout);
}

var exitWithError = (message) => {
  LOGGER.error( message );
  process.exit( 1 );
};

if ( token_config ) {
  let tokenContentConfig;
  try {
    let token_path = path.resolve(process.env.PWD, token_config);
    tokenContentConfig = JSON.parse(fs.readFileSync(token_path));
    if ( !tokenContentConfig.hasOwnProperty(TOKEN_CONFIG_KEY_NAME) ) {
      throw(new Error('cannot find token key in configuration'));
    }
    token = tokenContentConfig[TOKEN_CONFIG_KEY_NAME];
    LOGGER.info(`With token sourced from ${token_path}`);
  } catch ( e ) {
    exitWithError(`parsing configuration file ${e}`);
  }
  LOGGER.debug(`using token from configuration file`);
}

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
      // Data uris write to output and should not close it. We expect the output
      // to be closed only when the process finishes.
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
    if ( task.ca ) {
      options.ca = task.ca;
    }
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
        contentType = ( typeof contentType === 'string' ) ? contentType.toLowerCase()
                                                          : '';
        let parsedContentType = rangerRequest.parseContentType(contentType);
        if ( parsedContentType.json ) {
          if ( parsedContentType.version && !rangerRequest.supportedVersion(parsedContentType.version) ) {
            LOGGER.warn(
              `Unsupported streaming specification version in server response: ${parsedContentType.version}`
            );
          }
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
                output.end();
              };

              q.pause(); // To prevent run condition adding tasks vs processing queue
              /* jshint -W083 */
              // functions within a loop
              for ( var i = 0; i < uriData.uris.length; i++ ) {
                let newTask = {
                  uri:     uriData.uris[i],
                  headers: uriData.headers4uris[i]
                };
                if ( task.ca ) {
                  newTask.ca = task.ca;
                }
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
              q.resume();
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
              if ( trailerString ) {
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

          // Processing individual data uris should not try to close the output
          // stream. We expect the stream to be closed at the end of the whole
          // process.
          res.pipe(output, { end: false });
        }
      } else {
        let code = res.statusCode;
        let msg  = res.statusMessage || '';
        callback(`Non 200 status - ${code} ${msg}`);
      }
    });
  }
};

process.nextTick(() => {
  let task = { uri: url };
  if ( ca_file ) {
    LOGGER.debug(`Using ${ca_file} to try to source CAs.`);
    let ca_path = path.resolve(process.env.PWD, ca_file);
    let ca_content = fs.readFileSync(ca_path).toString();
    LOGGER.info(`With CA sourced from ${ca_path}`);
    task.ca = ca_content;
  }
  if ( token ) {
    let headers = {};
    headers[TOKEN_BEARER_KEY_NAME] = tokenUtils.formatTokenForHeader(token);
    task.headers = headers;
  }
  requestWorker(task, ( err ) => {
    if ( err ) {
      exitWithError( err );
    } else {
      LOGGER.debug('Success');
      process.exit(0);
    }
  });
});
