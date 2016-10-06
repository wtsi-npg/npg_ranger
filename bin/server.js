#!/usr/bin/env node

"use strict";

const fs      = require('fs');
const http    = require('http');
const assert  = require('assert');
const util    = require('util');
const MongoClient = require('mongodb').MongoClient;
const LOGGER      = require('../lib/logsetup.js');

const config = require('../lib/config.js');
// Call to config.provide() must occur here before requiring controller
// so that options object is built before it is provided to other modules.
const options = config.provide(config.fromCommandLine, true);
const RangerController = require('../lib/server/controller');

if ( options.get('debug') ) {
  LOGGER.level = 'debug';
}

LOGGER.info(config.logOpts());

/*
 * Main server script. Create the server object, establish database,
 * connection, setup server callbacks, start listening for incoming
 * requests.
 *
 * Providing config settings:
 *  Settings are provided from 3 locations:
 *  1. Command line - run with -h to see options.
 *  2. Config json file can be read if it is provided on command line
 *      by running with -c PATH or --configfile=PATH
 *  3. There are some defaults, which can be found in lib/config.js
 */

assert(process.env.USER, 'User environment variable is not defined');
const server = http.createServer();

// Exit gracefully on a signal to quit
process.on('SIGTERM', () => {
  server.close( () => {process.exit(0);} );
});
process.on('SIGINT', () => {
  server.close( () => {process.exit(0);} );
});

// Connect to the database and, if successful, define
// callbacks for the server.
var mongourl = options.get('mongourl');
MongoClient.connect(mongourl, options.get('mongoopt'), function(err, db) {

  assert.equal(err, null, `Failed to connect to ${mongourl}: ${err}`);
  LOGGER.info(`Connected to ${mongourl}`);

  var dbClose = (dbConn) => {
    if (dbConn) {
      LOGGER.info('Database connection closing');
      try {
        dbConn.close();
      } catch (err) {
        LOGGER.error(`Error closing db connection: ${err}`);
      }
    }
  };

  // Close database connection on server closing.
  server.on('close', () => {
    LOGGER.info("\nServer closing");
    dbClose(db);
  });

  // Exit gracefully on error, close the database
  // connection and remove the socket file.
  process.on('uncaughtException', (err) => {
    LOGGER.error(`Caught exception: ${err}\n`);
    dbClose(db);
    try {
      let port = options.get('port');
      if (typeof port != 'number') {
        // Throws an error if the assertion fails
        fs.accessSync(port, fs.W_OK);
        LOGGER.info(`Remove socket file ${port} that is left behind`);
        fs.unlinkSync(port);
      }
    } catch (err) {
      LOGGER.error(`Error removing socket file: ${err}`);
    }
    let code = 1;
    LOGGER.info(`Exiting with code ${code}`);
    process.exit(code);
  });

  // Set up a callback for requests.
  server.on('request', (request, response) => {
    if (options.get('debug')) {
      LOGGER.debug("MEMORY USAGE: " + util.inspect(process.memoryUsage()) + "\n");
    }

    // Ensure the processes initiated by request stops if the client disconnects.
    // Closing the response forces an error in the pipeline and allows for a
    // prompt closing of a socket established for this request.
    request.on('close', () => {
      LOGGER.info('CLIENT DISCONNECTED ');
      response.end();
    });

    // Create an instance of an application controller and let it
    // handle the request.
    LOGGER.debug('request headers: ' + JSON.stringify(request.headers));
    let controller = new RangerController(request, response, db);
    controller.handleRequest();
  });

  var createTempDataDir = () => {
    let tmpDir = options.get('tempdir');
    if (!fs.existsSync(tmpDir)) {
      fs.mkdirSync(tmpDir);
      LOGGER.debug(`Created temp data directory ${tmpDir}`);
    } else {
      LOGGER.debug(`Found temp data directory ${tmpDir}`);
    }
  };

  // Synchronously create directory for temporary data, then start listening.
  createTempDataDir(options.get('tempdir'));
  server.listen(options.get('port'), () => {
    LOGGER.info(`Server listening on ${options.get('hostname')}, ${options.get('port')}`);
  });
});
