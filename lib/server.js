"use strict";

const assert       = require('assert');
const EventEmitter = require('events');
const fs           = require('fs-extra');
const http         = require('http');
const MongoClient  = require('mongodb').MongoClient;
const util         = require('util');

const LOGGER           = require('../lib/logsetup.js');
const config           = require('../lib/config.js');
const RangerController = require('../lib/server/controller');

const SERVER_STARTED   = 'server_started';
const SERVER_CLOSED    = 'server_closed';

require('http-shutdown').extend();

class ServerFactory extends EventEmitter {

  /*
   * Server setup function.. Create the server object, establish database,
   * connection, setup server callbacks, start listening for incoming
   * requests.
   */
  startServer() {
    LOGGER.debug('Server factory start server');
    let options = config.provide();
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
        /*if (cluster.isWorker) {
          LOGGER.debug("WORKER: served by " + cluster.worker.id);
        }*/

        // Ensure the processes initiated by request stops if the client disconnects.
        // Closing the response forces an error in the pipeline and allows for a
        // prompt closing of a socket established for this request.
        request.on('close', () => {
          LOGGER.info('CLIENT DISCONNECTED');
          response.end();
        });

        // Create an instance of an application controller and let it
        // handle the request.
        let controller = new RangerController(request, response, db);
        controller.handleRequest(options.get('hostname'));
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
