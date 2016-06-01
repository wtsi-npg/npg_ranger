/* globals describe, it, expect */

"use strict";

const assert  = require('assert');
const os      = require('os');
const RangerModel = require('../lib/model.js');

describe('Class methods', function() {
  it('default format', function() {
    expect(RangerModel.defaultFormat()).toBe('bam');
  });
  it('supported formats', function() {
    expect(RangerModel.supportedFormats()).toEqual(['bam', 'cram', 'sam']);
  });
  it('is the format supported?', function() {
    expect( () => {RangerModel.supportsFormat();} )
      .toThrowError(assert.AssertionError,
      'Non-empty format string should be given');
    expect(RangerModel.supportsFormat('cram')).toBe(true);
    expect(RangerModel.supportsFormat('bed')).toBe(false);
  });
});

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
    expect( () => {m.process(1);} ).toThrowError(assert.AssertionError,
      'Query should be an object');
    expect( () => {m.process({one: 1});} ).toThrowError(assert.AssertionError,
      "Query should have a 'files' object inside");
    expect( () => {m.process({one: 1, files: 3});} )
      .toThrowError(assert.AssertionError,
      "Query should have a 'files' object inside");
    expect( () => {m.process({one: 1, files: {}});} )
      .toThrowError(assert.AssertionError,
      "Files object should have numeric 'length' property");
    expect( () => {m.process({one: 1, files: []});} )
      .toThrowError(assert.AssertionError,
      'Files should be given');
    expect( () => {m.process({files: ['file1', 'file2']});} ).toThrowError(assert.AssertionError,
      'Destination stream is required');
    expect( () => {m.process({files: ['file1', 'file2']}, {});} ).toThrowError(assert.AssertionError,
      'End callback is required');
    expect( () => {m.process({files: ['file1', 'file2']}, {}, {});} ).toThrowError(assert.AssertionError,
      'End callback is required');
  });  
});
