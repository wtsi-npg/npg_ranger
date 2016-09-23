/* globals describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, jasmine, spyOn */

"use strict";

const fse     = require('fs-extra');
const cluster = require('cluster');
const EventEmitter = require('events');
const config  = require('../../lib/config.js');
const assert  = require('assert');
// const RangerServer = require('../../lib/server.js');
const BinServer    = require('../../bin/server.js');

class EmptyServerFactory extends EventEmitter {
  startServer() {
    this.emit(BinServer.SERVER_STARTED);
  }
}

describe('Broker creation constructor validation', () => {
  it('constructor complained about missing serverFactory as param', () => {
    let bf = new BinServer.BrokerFactory();
    expect(() => { bf.buildBroker(); }).toThrowError(assert.AssertionError);
  });
});

describe('Decides which kind of broker', () => {
  let sf = new EmptyServerFactory();

  it('generates a flat broker when numworkers == 0', () => {
    let tmpDir = config.tempFilePath('npg_ranger_server_test_');
    config.provide(() => { return {
      tempdir: tmpDir,
      numworkers: 0
    };});

    let bf = new BinServer.BrokerFactory();
    let broker = bf.buildBroker(sf);
    expect(broker instanceof BinServer.Broker).toBe(true);
    expect(broker instanceof BinServer.FlatBroker).toBe(true);
    expect(broker instanceof BinServer.ClusteredBroker).toBe(false);
  });

  it('generates a cluster broker when requested one worker', () => {
    let tmpDir = config.tempFilePath('npg_ranger_server_test_');
    config.provide(() => { return {
      tempdir: tmpDir,
      numworkers: 1
    };});

    let bf = new BinServer.BrokerFactory();
    let broker = bf.buildBroker(sf);
    expect(broker instanceof BinServer.Broker).toBe(true);
    expect(broker instanceof BinServer.ClusteredBroker).toBe(true);
    expect(broker instanceof BinServer.FlatBroker).toBe(false);
  });

  it('Generates a cluster broker when requested two workers', () => {
    let tmpDir = config.tempFilePath('npg_ranger_server_test_');
    config.provide(() => { return {
      tempdir: tmpDir,
      numworkers: 2
    };});

    let bf = new BinServer.BrokerFactory();
    let broker = bf.buildBroker(sf);
    expect(broker instanceof BinServer.Broker).toBe(true);
    expect(broker instanceof BinServer.ClusteredBroker).toBe(true);
    expect(broker instanceof BinServer.FlatBroker).toBe(false);
  });
});

describe('Broker creation', () => {
  let brokerBuildCall;

  beforeAll( () => { });
  afterAll( () => { });

  beforeEach( () => {
    brokerBuildCall = jasmine.createSpy('brokerBuildCall');

    let factory = new BinServer.BrokerFactory();

    factory.on( BinServer.BROKER_BUILT, () => {
      brokerBuildCall();
    });

    factory.buildBroker(new EmptyServerFactory());
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

    let broker = new BinServer.BrokerFactory(new EmptyServerFactory()).buildBroker(new EmptyServerFactory());

    spyOn(cluster, 'fork');

    broker.on(BinServer.WORKER_FORKED, () => {
      workerForkedCall();
    });

    broker.on(BinServer.CLUSTER_STARTED, () => {
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

describe('Cluster limit consecutive forks', () => {
  // var sys = require('sys');
  var exec = require('child_process').exec;
  let child;

  afterEach(function() {
    if ( child.connected ) {
      child.disconnect();
    }
  });

  it('exits with correct code if max number of consec forks reached', ( done ) => {
    child = exec('bin/server.js -s -k 5 -l 1 -p 33000 -n10 -m mongodb://loclhost:27017/imc', (error) => {
      expect(error).not.toBe(null);
      expect(error.code).toEqual(210);
      done();
    });
  }, 10000);
});
