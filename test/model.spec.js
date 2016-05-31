/* globals describe, it, expect */

"use strict";

const assert  = require('assert');
const os      = require('os');
const RangerModel = require('../lib/model.js');

describe('Creating object instance', function() {
  it('temp directory attr is optional', function() {
    let m;
    expect( () => {m = new RangerModel();} ).not.toThrow();
    expect(m.tmpDir).toBe(os.tmpdir());
  });
  
  it('Temporary directory should exist', function() {
    expect( () => {new RangerModel('/some/dir');} )
      .toThrowError(assert.AssertionError,
      "Temp data directory '/some/dir' does not exist");
    let m;
    expect( () => {m = new RangerModel('test');} ).not.toThrow();
    expect(m.tmpDir).toBe('test');
  });
});

describe('Processing request', function() {
  let m = new RangerModel();
  it('Input validation', function() {
    expect( () => {m.process();} ).toThrowError(assert.AssertionError,
      'Query object is required');
    expect( () => {m.process({});} ).toThrowError(assert.AssertionError,
      'Destination stream is required');
    expect( () => {m.process({}, {});} ).toThrowError(assert.AssertionError,
      'End callback is required');
    expect( () => {m.process({}, {}, {});} ).toThrowError(assert.AssertionError,
      'End callback is required');
  });  
});
