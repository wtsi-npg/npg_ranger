"use strict";

const assert       = require('assert');
const EventEmitter = require('events');
const fs           = require('fs-extra');
const http         = require('http');
const MongoClient  = require('mongodb').MongoClient;
const util         = require('util');

const LOGGER           = require('winston');
const config           = require('../lib/config.js');
const RangerController = require('../lib/server/controller');

require('http-shutdown').extend();

/**
 * NodeJS implementation of util
 * @external util
 * @see  {@link https://nodejs.org/dist/latest-v4.x/docs/api/util.html|util}
 */

/**
  * fs-extra
  * @external fs-extra
  * @see {@link https://www.npmjs.com/package/fs-extra|fs-extra}
  */

/**
 * Official mongodb driver
 * @external mongodb
 * @see {@link https://www.npmjs.com/package/mongodb|mongodb}
 */

/**
 * Utility for gracefully close the server
 * @external http-shutdown
 * @see {@link https://www.npmjs.com/package/http-shutdown|http-shutdown}
 */

/**
 * <p>The business logic of the NPG Ranger server. The module is self sufficient
 * and can start a server itself in a simple format. A more powerful wrapper for
 * this module is provided at {@link module:bin/server|bin/server module}.</p>
 *
 * See {@link module:lib/server~ServerFactory|ServerFactory} for example use.
 * @module lib/server
 *
 * @requires {@link external:assert|assert}
 * @requires {@link external:events|events}
 * @requires {@link external:fs-extra|fs-extra}
 * @requires {@link external:http|http}
 * @requires {@link external:mongodb|mongodb}
 * @requires {@link external:util|util}
 * @requires {@link external:http-shutdown|http-shutdown}
 * @requires {@link external:winston|winston}
 *
 * @copyright Genome Research Limited 2016
 */

/**
 * Name of event emitted when server starts
 * @type {String}
 */
const SERVER_STARTED = 'server_started';
/**
 * Name of event emitted when server stops
 * @type {String}
 */
const SERVER_CLOSED  = 'server_closed';

/**
 *
 * @example
 *
 * let config       = require('../lib/config.js');
 * let RangerServer = require('./lib/server.js')
 *
 * // Assuming your server configuration comes from CLI options
 * config.provide(config.fromCommandLine);
 *
 * let serverFactory = new RangerServer.ServerFactory();
 * serverFactory.on(RangerServer.SERVER_STARTED, ( server ) => {
 *   console.log('server started');
 * });
 * serverFactory.on(RangerServer.SERVER_CLOSED, ( server ) => {
 *   console.log('server closed');
 * });
 * console.log('starting server');
 * serverFactory.startServer();
 *
 */
class ServerFactory extends EventEmitter {

  /**
   * Server setup function. Create the server object, establish database,
   * connection, setup server callbacks, start listening for incoming
   * requests.
   */
  startServer() {
    LOGGER.debug('Server factory start server');
    let options = config.provide();
    // Using withShutdown form http-shutdown
    let server = http.createServer().withShutdown();
    let self = this;

    // Exit gracefully on a signal to quit
    [ 'SIGTERM', 'SIGINT', 'SIGHUP' ].forEach( ( sig ) => {
      process.on( sig, () => {
        server.shutdown( () => {
          process.exit(0);
        });
      });
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
        self.emit(SERVER_CLOSED, server);
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
        LOGGER.debug("MEMORY USAGE: " + util.inspect(process.memoryUsage()) + "\n");

        // Ensure the processes initiated by request stops if the client disconnects.
        // Closing the response forces an error in the pipeline and allows for a
        // prompt closing of a socket established for this request.
        request.on('close', () => {
          LOGGER.info('CLIENT DISCONNECTED');
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
        LOGGER.debug(`Using temp data directory ${tmpDir}`);
        fs.ensureDirSync(tmpDir);
      };

      // Synchronously create directory for temporary data, then start listening.
      createTempDataDir(options.get('tempdir'));
      server.listen(options.get('port'), () => {
        LOGGER.info(`Server listening on ${options.get('hostname')}, ${options.get('port')}`);
        self.emit(SERVER_STARTED, server);
      });
    });
  }
}

module.exports = {
  ServerFactory:  ServerFactory,
  SERVER_STARTED: SERVER_STARTED,
  SERVER_CLOSED:  SERVER_CLOSED
};
