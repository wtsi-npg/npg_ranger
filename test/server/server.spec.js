/* globals describe, it, expect, beforeAll, afterAll, beforeEach, jasmine, spyOn */

"use strict";

const fse     = require('fs-extra');
const cluster = require('cluster');
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

    let broker = new Server.BrokerFactory().buildBroker(numworkers);

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
