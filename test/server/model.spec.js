/* globals describe, it, expect, beforeAll */

"use strict";

const assert  = require('assert');
const fs      = require('fs');
const RangerModel = require('../../lib/server/model.js');
const config = require('../../lib/config.js');

// Create temp dir here so it is available for all tests.
// Use this dir as a default dir that will be available in all.
var tmpDir    = config.tempFilePath('npg_ranger_model_test_');
if (!fs.existsSync(tmpDir)) {
  fs.mkdirSync(tmpDir);
}
var dummy     = function() { return {tempdir: tmpDir}; };
config.provide(dummy);

describe('Class methods', function() {
  beforeAll(function() {
    config.provide(dummy);
  });
  it('default format', function() {
    expect(RangerModel.defaultFormat()).toBe('BAM');
  });
  it('supported formats', function() {
    expect(RangerModel.supportedFormats()).toEqual(['BAM', 'CRAM', 'SAM']);
  });
  it('is the format supported?', function() {
    expect( () => {RangerModel.supportsFormat();} )
      .toThrowError(assert.AssertionError,
      'Non-empty format string should be given');
    expect(RangerModel.supportsFormat('CRAM')).toBe(true);
    expect(RangerModel.supportsFormat('bed')).toBe(false);
    expect(RangerModel.supportsFormat('BED')).toBe(false);
  });
  it('textual formats', function() {
    expect(RangerModel.textualFormats()).toEqual(['SAM']);
  });
  it('is the format textual?', function() {
    expect( () => {RangerModel.isTextualFormat();} )
      .toThrowError(assert.AssertionError,
      'Non-empty format string should be given');
    expect(RangerModel.isTextualFormat('CRAM')).toBe(false);
    expect(RangerModel.isTextualFormat('BED')).toBe(false);
    expect(RangerModel.isTextualFormat('SAM')).toBe(true);
  });
});

describe('Creating object instance', function() {
  beforeAll(function() {
    config.provide(dummy);
  });
  it('temp directory attr is optional', function() {
    let m;
    expect( () => {m = new RangerModel();} ).not.toThrow();
    expect(m.tmpDir).toBe(tmpDir);
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
   beforeAll(function() {
    config.provide(dummy);
  });
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
