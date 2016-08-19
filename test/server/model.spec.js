/* globals describe, it, expect */

"use strict";

const assert  = require('assert');
const os      = require('os');
const RangerModel = require('../../lib/server/model.js');

describe('Class methods', function() {
  it('default format', function() {
    expect(RangerModel.defaultFormat()).toBe('BAM');
  });
  it('supported formats', function() {
    expect(RangerModel.supportedFormats()).toEqual(['BAM', 'CRAM', 'SAM', 'VCF']);
  });
  it('is the format supported?', function() {
    expect( () => {RangerModel.supportsFormat();} )
      .toThrowError(assert.AssertionError,
      'Non-empty format string should be given');
    expect(RangerModel.supportsFormat('CRAM')).toBe(true);
    expect(RangerModel.supportsFormat('VCF')).toBe(true);
    expect(RangerModel.supportsFormat('bed')).toBe(false);
    expect(RangerModel.supportsFormat('BED')).toBe(false);
  });
  it('textual formats', function() {
    expect(RangerModel.textualFormats()).toEqual(['SAM', 'VCF']);
  });
  it('is the format textual?', function() {
    expect( () => {RangerModel.isTextualFormat();} )
      .toThrowError(assert.AssertionError,
      'Non-empty format string should be given');
    expect(RangerModel.isTextualFormat('CRAM')).toBe(false);
    expect(RangerModel.isTextualFormat('BED')).toBe(false);
    expect(RangerModel.isTextualFormat('VCF')).toBe(true);
    expect(RangerModel.isTextualFormat('SAM')).toBe(true);
  });
});

describe('Creating object instance', function() {
  it('temp directory attr is optional', function() {
    let m;
    expect( () => {m = new RangerModel(1000);} ).not.toThrow();
    expect(m.tmpDir).toBe(os.tmpdir());
  });

  it('Temporary directory should exist', function() {
    expect( () => {new RangerModel(2000, '/some/dir');} )
      .toThrowError(assert.AssertionError,
      "Temp data directory '/some/dir' does not exist");
    let m;
    expect( () => {m = new RangerModel(1000, 'test');} ).not.toThrow();
    expect(m.tmpDir).toBe('test');
  });
});

describe('Processing request', function() {
  let m = new RangerModel(1000);
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
    expect( () => {m.process({files: ['file1'], format: 'VCF'}, {}, () => {return;});} ).toThrowError(ReferenceError,
      'database does not hold location of reference .fa file');
  });
});
