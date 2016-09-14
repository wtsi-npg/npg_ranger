#!/usr/bin/env node

"use strict";

const fs      = require('fs-extra');
const http    = require('http');
const assert  = require('assert');
const util    = require('util');
const EventEmitter = require('events');
const MongoClient  = require('mongodb').MongoClient;
const LOGGER       = require('../lib/logsetup.js');

const config = require('../lib/config.js');

class RangerBroker extends EventEmitter {
  startServer( options ) {
    const server = http.createServer();
    let broker = this;

    // Exit gracefully on a signal to quit
    [ 'SIGTERM', 'SIGINT', 'SIGHUP' ].forEach( ( sig ) => {
      process.on( sig, () => {
        server.close( () => {
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
        broker.emit('workerClose');
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

        const RangerController = require('../lib/server/controller');

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
        broker.emit('serverStarted', server);
      });
    });
  }
}

class FlatBroker extends RangerBroker {
  start( options ) {
    this.startServer(options);
  }
}

class ClusteredBroker extends RangerBroker {
  start( options ) {
    const cluster = require('cluster');
    if ( cluster.isMaster ) {
      LOGGER.info(config.logOpts());
      for (let i = 0; i < options.get('numworkers'); i++) {
        this.emit('workerForked');
        cluster.fork();
      }
      let consec = 0;
      cluster.on('exit', (worker, code, signal) => {
        LOGGER.debug('Worker %d died (%s). Forking to replace ...', worker.id, signal || code);
        let waitingConsec = options.get('clustertimeout');
        let maxConsec = options.get('clustermaxdeaths');
        LOGGER.debug(`${consec} forks have died in the previous ${waitingConsec} seconds.`);
        if ( consec >= maxConsec ) {
          LOGGER.error('Too many forks started in short span of time. Trying to exit now.');
          process.exit(1);
        }
        consec += 1;
        cluster.fork();
        this.emit('workerForked');
        setTimeout( () => {
          consec -= 1;
        }, waitingConsec * 1000 );
      });
      this.emit('clusterStarted', cluster);
    } else {
      if ( cluster.isWorker ) {
        this.emit('workerStarted', cluster.worker);
        LOGGER.debug('WORKER: new fork ' + cluster.worker.id);
        this.startServer(options);
      }
    }
  }
}

class BrokerFactory extends EventEmitter {
  buildBroker(numWorkers) {
    // assert( typeof numWorkers === 'number', 'Number of workers must be a number' );
    let broker;
    if ( numWorkers === 0 ) {
      broker = new FlatBroker();
    } else {
      broker = new ClusteredBroker();
    }
    this.emit('brokerBuilt', broker);
    return broker;
  }
}

if ( require.main === module ) {
  // Call to config.provide() must occur here before requiring controller
  // so that options object is built before it is provided to the
  // other modules.
  const options = config.provide(config.fromCommandLine);

  const numWorkers = options.get('numworkers');

  if ( options.get('debug') ) {
    LOGGER.level = 'debug';
  }

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
  let bf = new BrokerFactory();
  let broker = bf.buildBroker(numWorkers);

  broker.on('clusterStarted', () => { console.log('cluster started'); });
  broker.on('workerStarted', () => { console.log('worker started'); });
  broker.on('workerForked', () => { console.log('worker forked'); });

  broker.start(options);
}
