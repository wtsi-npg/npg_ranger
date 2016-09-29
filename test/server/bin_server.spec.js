/* globals describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, jasmine, spyOn, fail */

"use strict";

const fse          = require('fs-extra');
const cluster      = require('cluster');
const EventEmitter = require('events');
const assert       = require('assert');
const tmp          = require('tmp');

const config       = require('../../lib/config.js');
const BinServer = require('../../bin/server.js');

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

describe('Flat server', () => {
  let server;
  let http = require('http');
  let tmpDir = config.tempFilePath('npg_ranger_server_test_');
  let testConfBuilder = () => {
    return {
      tempdir: tmpDir,
      numworkers: 0
    };
  };
  let content = 'some chars';

  class SimpleServerFactory extends EventEmitter {
    startServer() {
      server = http.createServer( (request, response) => {
        response.writeHead(200, {"Content-Type": "text/html"});
        response.write(content);
        response.end();
      });
      server.listen(0, () => {
        this.emit(BinServer.SERVER_STARTED);
      });
    }
  }

  beforeAll( () => {
    fse.ensureDirSync(tmpDir);
    config.provide(testConfBuilder);
  });

  afterAll( () => {
    fse.removeSync(tmpDir);
  });

  beforeEach( (done) => {
    let factory = new SimpleServerFactory();
    let broker = new BinServer.BrokerFactory().buildBroker(factory);
    factory.on(BinServer.SERVER_STARTED, () => {
      done();
    });
    broker.start();
  });

  afterEach(( done ) => {
    server.close(() => {
      done();
    });
  });

  it('starts a test http server', ( done ) => {
    let requestURL = 'http://localhost:' + server.address().port + '/';
    http.get(requestURL, ( res ) => {
      expect(res.statusCode).toBe(200);
      let responseContent = '';
      res.on('data', (data) => {
        responseContent += data;
      });
      res.on('end', () => {
        expect(responseContent).toBe(content);
        done();
      });
    });
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

  beforeAll( () => {
    fse.ensureDirSync(tmpDir);
    config.provide(testConfBuilder);
    oldFork = cluster.fork;
  });

  afterAll( () => {
    fse.removeSync(tmpDir);
    cluster.fork = oldFork;
  });

  beforeEach( () => {
    clusterStartedCall = jasmine.createSpy('clusterStartedCall');
    workerForkedCall = jasmine.createSpy('workerForkedCall');

    let broker = new BinServer.BrokerFactory().buildBroker(new EmptyServerFactory());

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

describe('Cluster limit consecutive forks at start', () => {
  let exec = require('child_process').exec;
  let child;

  afterEach( () => {
    if ( child.connected ) {
      child.disconnect();
    }
  });

  it('exits with correct code if max number of consec forks reached', ( done ) => {
    child = exec('bin/server.js -k5 -l1 -n10 -m mongodb://loclhost:27017/imc', (error) => {
      expect(error).not.toBe(null);
      if ( !!error ) {
        expect(error.code).toEqual(210);
      }
      done();
    });
  }, 10000);
});

describe('Cluster limit consecutive forks', () => {
  let spawn    = require('child_process').spawn;
  let exec     = require('child_process').exec;
  let execSync = require('child_process').execSync;
  let child;

  let tmp_dir;

  let BASE_PORT  = 1400;
  let PORT_RANGE = 200;
  let PORT = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;

  beforeAll(() => {
    let tmpobj = tmp.dirSync({ prefix: 'npg_ranger_test_' });
    tmp_dir = tmpobj.name;
    let command = `mongod -f test/server/data/mongodb_conf.yml --port ${PORT} --dbpath ${tmp_dir} --pidfilepath ${tmp_dir}/mpid --logpath ${tmp_dir}/dbserver.log`;
    let out = execSync(command);
    console.log(`Loaded data to MONGO DB: ${out}`);
  });

  afterAll(() => {
    execSync(`mongo 'mongodb://localhost:${PORT}/admin' --eval 'db.shutdownServer()'`);
    fse.remove(tmp_dir, (err) => {if (err) {console.log(`Error removing ${tmp_dir}: ${err}`);}});
  });

  afterEach( () => {
    child.kill('SIGINT');
  });

  it('re forks if single child dies', (done) => {
    let command = 'bin/server.js';
    let numForks = 7;
    child = spawn(command,[
      `-mmongodb://localhost:${PORT}`,
      '-k5', '-l1', // Max 3 deaths in 1 second
      '-n' + numForks,
      '-p33000']
    );
    setTimeout(() => {
      let grandchildrenBefore, grandchildrenAfter;
      exec('pgrep -P ' + child.pid, (error, stdout) => {
        stdout = stdout.trim();
        grandchildrenBefore = stdout.split('\n').map((value) => {
          return Number.parseInt(value);
        });
        let victim = grandchildrenBefore[0];
        expect(grandchildrenBefore.length).toEqual(numForks);
        process.kill(victim, 'SIGKILL');
        setTimeout(() => {
          exec('pgrep -P ' + child.pid, (error2, stdout2) => {
            stdout2 = stdout2.trim();
            grandchildrenAfter = stdout2.split('\n').map((value) => {
              return Number.parseInt(value);
            });
            expect(grandchildrenAfter.length).toEqual(numForks);
            grandchildrenBefore.slice(1).forEach((value)=>{
              expect(grandchildrenAfter.indexOf(value)).not.toEqual(-1);
            });
            expect(grandchildrenAfter.indexOf(victim)).toEqual(-1);
            done();
          });
        }, 1000);
      });
    }, 3000);
  }, 7000);

  it('dies with correct error code if enough children die in short time', ( done ) => {
    let command = 'bin/server.js';
    let killing = false;
    child = spawn(command,[
      `-mmongodb://localhost:${PORT}`,
      '-k3', '-l1', // Max 3 deaths in 1 second
      '-n5',
      '-p33000']
    );
    child.on('close', (code) => {
      expect(killing).toBe(true);
      expect(code).toEqual(210);
      done();
    });
    let grandchildren;
    setTimeout(() => {
      exec('pgrep -P ' + child.pid, (error, stdout) => {
        killing = true;
        stdout = stdout.trim();
        grandchildren = stdout.split('\n').map((value) => {
          return Number.parseInt(value);
        });
        grandchildren.forEach((value) => {
          try {
            process.kill(value, 'SIGKILL');
          } catch(e) {
            fail('Something went wrong when killing grandchildren');
          }
        });
      });
    }, 5000);
  }, 10000);
});
