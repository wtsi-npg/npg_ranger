/* globals describe, it, expect, beforeEach */

"use strict";

const spawn    = require('child_process').spawn;
const devnull  = require('dev-null');
const pipeline = require('../../lib/server/pipeline.js');

var path = 'test/server/data/pipeline/text.txt';
var isSuccess = null;

var setup = function(done, cat_command, f, wc_command, wc_options) {

  const cat = spawn(cat_command, [f]);
  cat.title = 'cat';
  var prs = [cat];

  if (wc_command) {
    if (! wc_options) {
      wc_options = '-l';
    }
    const wc  = spawn(wc_command, [wc_options]);
    wc.title = 'wc';
    prs.push(wc);
  }

  // The following two functions will be called asynchronously
  isSuccess = null;
  var success = function() {isSuccess = true;  done();};
  var failure = function() {isSuccess = false; done();};

  var pline = pipeline(prs, success, failure);
  pline.run(devnull());
};

describe('Run function input validation', function() {

  const cat = spawn('cat', [path]);
  var fun = function() {};

  it('Destination is not given - error', function() {
    expect( () => {pipeline([cat], fun, fun).run();} ).toThrowError(ReferenceError,
    'Destination stream is not defined');
  });


  it('Array of processes should be defined', function() {
    expect( () => {pipeline().run(devnull());} ).toThrowError(
      ReferenceError, 'Array of processes should be defined'
    );
  });
  it('Array of processes cannot be null', function() {
    expect( () => {pipeline(null, fun, fun).run(devnull());} ).toThrowError(
      ReferenceError, 'Array of processes should be defined'
    );
  });
  it('Processes should be an array', function() {
    expect( () => {pipeline(fun, fun).run(devnull());} ).toThrowError(
      TypeError, 'processes should be an array'
    );
  });
  it('Processes array cannot be empty', function() {
    expect( () => {pipeline([], fun, fun).run(devnull());} ).toThrowError(
      RangeError, 'processes array cannot be empty'
    );
  });
  it('All process array members have to be defined', function() {
    expect( () => {pipeline([cat, null, cat], fun, fun).run(devnull());} ).toThrowError(
      ReferenceError, 'Undefined process at index 1'
    );
  });
  it('Process shoudl be an instance of EventEmitter', function() {
    expect( () => {pipeline([cat, 'dog', cat], fun, fun).run(devnull());} ).toThrowError(
      TypeError, 'Not an event emitter at index 1'
    );
  });

  it('Success callback should be defined', function() {
    expect( () => {pipeline([cat]).run(devnull());} ).toThrowError(
      ReferenceError, 'Success callback should be defined'
    );
  });
  it('Success callback should be a function', function() {
    expect( () => {pipeline([cat], 'callback').run(devnull());} ).toThrowError(
      TypeError, 'Success callback should be a function'
    );
  });
  it('Failure callback should be defined', function() {
    expect( () => {pipeline([cat], fun).run(devnull());} ).toThrowError(
      ReferenceError, 'Failure callback should be defined'
    );
  });
  it('Failure callback should be a function', function() {
    expect( () => {pipeline([cat], fun, 'callback').run(devnull());} ).toThrowError(
      TypeError, 'Failure callback should be a function'
    );
  });
});

describe('Valid two-process pipeline', function() {
  beforeEach(function(done) {
    setup(done, 'cat', path, 'wc');
  });
  it('Pipeline succeeds', function(done) {
    expect(isSuccess).toBe(true);
    done();
  });
});

describe('Valid one-process pipeline', function() {
  beforeEach(function(done) {
    setup(done, 'cat', path);
  });
  it('Pipeline succeeds', function(done) {
    expect(isSuccess).toBe(true);
    done();
  });
});

describe('Valid two-process pipeline, invalid input', function() {
  beforeEach(function(done) {
    setup(done, 'cat', 'test/data/pipeline/test.txt', 'wc');
  });
  it('Pipeline fails', function(done) {
    expect(isSuccess).toBe(false);
    done();
  });
});

describe('Valid one-process pipeline, invalid input', function() {
  beforeEach(function(done) {
    setup(done, 'cat', 'test/data/pipeline/test.txt');
  });
  it('Pipeline fails', function(done) {
    expect(isSuccess).toBe(false);
    done();
  });
});

describe('Invalid pipeline, first command incorrect', function() {
  beforeEach(function(done) {
    setup(done, 'notcat', path, 'wc');
  });
  it('Pipeline fails', function(done) {
    expect(isSuccess).toBe(false);
    done();
  });
});

describe('Invalid pipeline, second command incorrect', function() {
  beforeEach(function(done) {
    setup(done, 'cat', path, 'nowc');
  });
  it('Pipeline fails', function(done) {
    expect(isSuccess).toBe(false);
    done();
  });
});

describe('Invalid pipeline, second command option is invalid', function() {
  beforeEach(function(done) {
    setup(done, 'cat', path, 'wc', '-something');
  });
  it('Pipeline fails', function(done) {
    expect(isSuccess).toBe(false);
    done();
  });
});
