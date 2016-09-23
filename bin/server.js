#!/usr/bin/env node

"use strict";

const EventEmitter = require('events');
const assert       = require('assert');
const LOGGER       = require('../lib/logsetup.js');

const config       = require('../lib/config.js');
const RangerServer = require('../lib/server.js');

const SERVER_STARTED  = RangerServer.SERVER_STARTED;
const SERVER_CLOSED   = RangerServer.SERVER_CLOSED;

const BROKER_BUILT          = 'broker_built';
const CLUSTER_STARTED       = 'cluster_started';
const WORKER_STARTED        = 'worker_started';
const WORKER_FORKED         = 'worker_forked';
const WORKER_CLOSED         = 'worker_closed';
const HIT_LIMIT_CONSEC_FORK = 'hit_limit_consec_forks';

const ERROR_SERVER_LIMIT_CONSEC_FORK = 210;

class Broker extends EventEmitter {
  constructor(serverFactory) {
    super();
    this.serverFactory = serverFactory;
    this.serverFactory.on(SERVER_STARTED, ( server ) => {
      this.emit(SERVER_STARTED, server);
    });
    this.serverFactory.on(SERVER_CLOSED, ( server ) => {
      this.emit(SERVER_CLOSED, server);
    });
  }

  start() {
    LOGGER.debug('broker start server');
  }
}

class FlatBroker extends Broker {
  start() {
    super.start();
    this.serverFactory.startServer();
  }
}


class ClusteredBroker extends Broker {
  start() {
    super.start();
    let options = config.provide();
    const cluster = require('cluster');
    if ( cluster.isMaster ) {
      LOGGER.info(config.logOpts());
      for (let i = 0; i < options.get('numworkers'); i++) {
        cluster.fork();
        this.emit(WORKER_FORKED);
      }
      let consec = 0;
      let exiting = false;
      cluster.on('exit', (worker, code, signal) => {
        LOGGER.debug('Worker %d died (%s). Forking to replace ...', worker.id, signal || code);
        let waitingConsec = options.get('clustertimeout');
        let maxConsec = options.get('clustermaxdeaths');
        LOGGER.debug(`${consec} forks have died in the previous ${waitingConsec} seconds.`);
        if ( consec >= maxConsec ) {
          if ( !exiting ) {
            exiting = true;
            LOGGER.error('Too many forks started in short span of time. Trying to exit now.');
            this.emit(HIT_LIMIT_CONSEC_FORK);
            setTimeout( () => {
              process.exit(ERROR_SERVER_LIMIT_CONSEC_FORK);
            }, 3000);
          }
        } else {
          consec += 1;
          cluster.fork();
          this.emit(WORKER_FORKED);
          setTimeout( () => {
            consec -= 1;
          }, waitingConsec * 1000 );
        }
      });
      this.emit(CLUSTER_STARTED, cluster);
    } else {
      if ( cluster.isWorker ) {
        this.emit(WORKER_STARTED, cluster.worker);
        LOGGER.debug('WORKER: new fork ' + cluster.worker.id);
        this.serverFactory.startServer();
      }
    }
  }
}

class BrokerFactory extends EventEmitter {
  buildBroker(serverFactory) {
    assert(serverFactory, 'serverFactory is required');
    let options = config.provide();
    let numworkers = Number.parseInt(options.get('numworkers'));
    assert(Number.isInteger(numworkers), 'numworkers must be an integer number');
    let broker;
    if ( !numworkers ) {
      broker = new FlatBroker(serverFactory);
    } else {
      broker = new ClusteredBroker(serverFactory);
    }
    this.emit(BROKER_BUILT, broker);
    return broker;
  }
}

if ( require.main === module ) {
  // Providing config settings:
  //  Settings are provided from 3 locations:
  //  1. Command line - run with -h to see options.
  //  2. Config json file can be read if it is provided on command line
  //      by running with -c PATH or --configfile=PATH
  //  3. There are some defaults, which can be found in lib/config.js
  // Call to config.provide() must occur here before requiring controller
  // so that options object is built before it is provided to the
  // other modules
  const options = config.provide(config.fromCommandLine);

  if ( options.get('debug') ) {
    LOGGER.level = 'debug';
  }

  assert(process.env.USER, 'User environment variable is not defined');
  let bf = new BrokerFactory();
  let sf = new RangerServer.ServerFactory();

  sf.on(SERVER_STARTED, () => { LOGGER.debug('require.main sf.server_started'); });
  sf.on(SERVER_CLOSED,  () => { LOGGER.debug('require.main sf.server_closed'); });

  let broker = bf.buildBroker(sf);

  broker.on(CLUSTER_STARTED, () => { LOGGER.debug('require.main cluster started'); });
  broker.on(WORKER_STARTED,  () => { LOGGER.debug('require.main worker started'); });
  broker.on(WORKER_FORKED,   () => { LOGGER.debug('require.main worker forked'); });

  process.nextTick(() => {
    broker.start();
  });
}

module.exports = {
  BrokerFactory:   BrokerFactory,
  Broker:          Broker,
  FlatBroker:      FlatBroker,
  ClusteredBroker: ClusteredBroker,
  BROKER_BUILT:    BROKER_BUILT,
  CLUSTER_STARTED: CLUSTER_STARTED,
  WORKER_STARTED:  WORKER_STARTED,
  WORKER_FORKED:   WORKER_FORKED,
  WORKER_CLOSED:   WORKER_CLOSED,
  SERVER_STARTED:  SERVER_STARTED,
  SERVER_CLOSED:   SERVER_CLOSED
};
