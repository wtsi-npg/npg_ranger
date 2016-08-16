#!/usr/bin/env node

"use strict";

const os      = require('os');
const fs      = require('fs');
const path    = require('path');
const http    = require('http');
const assert  = require('assert');
const util    = require('util');
const MongoClient = require('mongodb').MongoClient;
const GetOpt      = require('node-getopt');
const LOGGER      = require('../lib/logsetup.js');

const RangerController = require('../lib/server/controller');

var opt = new GetOpt([
    ['p','port=PORT'        ,'PORT or socket which server listens on'],
    ['m','mongourl=URI'     ,'URI to contact mongodb'],
    ['t','tempdir=PATH'     ,'PATH of temporary directory'],
    ['H','hostname=HOST'    ,'override hostname with HOST'],
    ['S','strictmode'       ,'run in strict mode'],
    ['s','skipauth'         ,'skip authorisation steps'],
    ['d','debug'            ,'debugging mode for this server'],
    ['h','help'             ,'display this help']
]).bindHelp().parseSystem();

const PORT               = opt.options.port || opt.argv[0]
                           || path.join(os.tmpdir(), process.env.USER, 'npg_ranger.sock');
const HOST               = opt.options.hostname || os.hostname() || 'localhost';
const MONGO              = opt.options.mongourl || 'mongodb://sf2-farm-srv1:27017/imetacache';
const TEMP_DATA_DIR_NAME = 'npg_ranger_data';
const TEMP_DATA_DIR      = opt.options.tempdir || path.join(os.tmpdir(), process.env.USER, TEMP_DATA_DIR_NAME);

const MONGO_OPTIONS = {
  db: {
    numberOfRetries: 5
  },
  server: {
    auto_reconnect: true,
    poolSize: 40,
    socketOptions: {
      connectTimeoutMS: 5000
    }
  },
  replSet: {},
  mongos: {}
};

if ( opt.options.debug ) {
  LOGGER.level = 'debug';
}
LOGGER.info(opt.options);

/*
 * Main server script. Create the server object, establish database,
 * connection, setup server callbacks, start listening for incoming
 * requests.
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
MongoClient.connect(MONGO, MONGO_OPTIONS, function(err, db) {

  assert.equal(err, null, `Failed to connect to ${MONGO}: ${err}`);
  LOGGER.info(`Connected to ${MONGO}`);

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
      if (typeof PORT != 'number') {
        // Throws an error if the assertion fails
        fs.accessSync(PORT, fs.W_OK);
        LOGGER.info(`Remove socket file ${PORT} that is left behind`);
        fs.unlinkSync(PORT);
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
    if (opt.options.debug) {
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
    let controller = new RangerController(
      request, response, db, TEMP_DATA_DIR, opt.options.skipauth,
      opt.options.strict);
    controller.handleRequest(HOST);
  });

  var createTempDataDir = (tmpDir) => {
    if (!fs.existsSync(tmpDir)) {
      let dir = path.join(os.tmpdir(), process.env.USER);
      if (!fs.existsSync(dir)) {
        fs.mkdirSync(dir);
      }
      fs.mkdirSync(tmpDir);
      LOGGER.debug(`Created temp data directory ${tmpDir}`);
    } else {
      LOGGER.debug(`Found temp data directory ${tmpDir}`);
    }
  };

  // Synchronously create directory for temporary data, then start listening.
  createTempDataDir(TEMP_DATA_DIR);
  server.listen(PORT, () => {
    LOGGER.info(`Server listening on ${HOST}, ${PORT}`);
  });
});

