/* globals describe, it, expect, beforeAll, afterAll, beforeEach, jasmine, spyOn */

"use strict";

const fse     = require('fs-extra');
const cluster = require('cluster');
const EventEmitter = require('events');
const config  = require('../../lib/config.js');
const Server  = require('../../bin/server.js');

describe('Broker creation', () => {
  let brokerBuildCall;

  beforeAll( () => { });
  afterAll( () => { });

  beforeEach( () => {
    brokerBuildCall = jasmine.createSpy('brokerBuildCall');

    let factory = new Server.BrokerFactory();

    factory.on( Server.BROKER_BUILT, () => {
      brokerBuildCall();
    });

    factory.buildBroker();
  });

  it('emitted the build broker event', ( done ) => {
    expect(brokerBuildCall).toHaveBeenCalled();
    done();
  });

  it('emitted the build broker event once', ( done ) => {
    expect(brokerBuildCall.calls.count()).toEqual(1);
    done();
  });
});

describe('Cluster creation', () => {
  let clusterStartedCall;
  let workerForkedCall;
  let numworkers = 3;
  let oldFork;

  let tmpDir = config.tempFilePath('npg_ranger_server_test_');
  let testConfBuilder = () => {
    return {
      tempdir: tmpDir,
      numworkers: numworkers
    };
  };

  beforeAll(function() {
    fse.ensureDirSync(tmpDir);
    config.provide(testConfBuilder);
    oldFork = cluster.fork;
  });

  afterAll(function() {
    fse.removeSync(tmpDir);
    cluster.fork = oldFork;
  });

  beforeEach( () => {
    clusterStartedCall = jasmine.createSpy('clusterStartedCall');
    workerForkedCall = jasmine.createSpy('workerForkedCall');

    let broker = new Server.BrokerFactory().buildBroker();

    spyOn(cluster, 'fork');

    broker.on(Server.WORKER_FORKED, () => {
      workerForkedCall();
    });

    broker.on(Server.CLUSTER_STARTED, () => {
      clusterStartedCall();
    });

    broker.start();
  });

  it('emitted the forked worker event', ( done ) => {
    expect(workerForkedCall).toHaveBeenCalled();
    done();
  });

  it('emitted the forked worker event correct number ot times', ( done ) => {
    expect(workerForkedCall.calls.count()).toEqual(numworkers);
    done();
  });

  it('caused fork to be called correct number ot times', ( done ) => {
    expect(cluster.fork.calls.count()).toEqual(numworkers);
    done();
  });

  it('emitted the cluster start event', ( done ) => {
    expect(clusterStartedCall).toHaveBeenCalled();
    done();
  });

  it('emitted the cluster start event once', ( done ) => {
    expect(clusterStartedCall.calls.count()).toEqual(1);
    done();
  });
});

describe('Other test', () => {
  let tmpDir = config.tempFilePath('npg_ranger_server_test_');
  let numworkers = 3;
  let oldFork;
  let testConfBuilder = () => {
    return {
      tempdir: tmpDir,
      numworkers: numworkers
    };
  };
  let workers;

  var FakeWorker = class FakeWorker extends EventEmitter {
    constructor(id) {
      super();
      this.id = id;
      this.dead = false;
    }

    justDie(exitCode, signalCode) {
      delete workers['' + this.id];
      this.suicide = true;
      this.state = 'dead';
      this.emit('exit', exitCode, signalCode);
      cluster.emit('exit', this, exitCode, signalCode);
    }

    disconnect() {
      this.dead = true;
      this.suicide = true;
    }

    kill() {
      this.dead = true;
      this.suicide = false;
    }

    isDead() {
      return this.dead;
    }
  };

  beforeAll(function() {
    fse.ensureDirSync(tmpDir);
    config.provide(testConfBuilder);
    oldFork = cluster.fork;
  });

  afterAll(function() {
    fse.removeSync(tmpDir);
    cluster.fork = oldFork;
  });

  beforeEach( () => {
    let broker = new Server.BrokerFactory().buildBroker();
    let ids = 0;
    workers = {};
    cluster.workers = {};

    spyOn(cluster, 'fork').and.callFake(function() {
      cluster.setupMaster();
      let fakeWorker = new FakeWorker(ids++);
      console.log('Fork ' + ids);

      workers[fakeWorker.id] = fakeWorker;
      cluster.workers[fakeWorker.id] = fakeWorker;

      fakeWorker.on('message', (message, handle) =>
        cluster.emit('message', message, handle)
      );

      // let emitForkNT = function(worker) {
      cluster.emit('fork', fakeWorker);
      // };
      // process.nextTick(emitForkNT, fakeWorker);
      return fakeWorker;
    });

    broker.start();
  });

  it('caused fork to be called correct number ot times', ( done ) => {
    expect(cluster.fork.calls.count()).toEqual(numworkers);
    // workers['0'].justDie(1, null);
    expect(cluster.fork.calls.count()).toEqual(numworkers);
    done();
  });
});
