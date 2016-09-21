/* globals describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, jasmine, spyOn */

"use strict";

const fse     = require('fs-extra');
const cluster = require('cluster');
// const EventEmitter = require('events');
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
    child = exec('bin/server.js -s -k 5 -l 1 -p 33000 -n4 -m mongodb://loclhost:27017/imc', (error) => {
      expect(error).not.toBe(null);
      expect(error.code).toEqual(210);
      done();
    });
  }, 10000);
});