/* globals describe, it, expect, beforeAll, afterAll */

"use strict";

const assert  = require('assert');
const fse     = require('fs-extra');
const RangerModel = require('../../lib/server/model.js');
const config = require('../../lib/config.js');

// Create temp dir here so it is available for all tests.
// Use this dir as a default dir that will be available in all.
var tmpDir    = config.tempFilePath('npg_ranger_model_test_');
var dummy     = function() { return {tempdir: tmpDir}; };
config.provide(dummy);

describe('Class methods', function() {
  beforeAll(function() {
    fse.ensureDirSync(tmpDir);
    config.provide(dummy);
  });

  afterAll(function() {
    fse.removeSync(tmpDir);
  });

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
  it('does the query have a reference?', () => {
    expect( () => {RangerModel.hasReference();} )
      .toThrowError(assert.AssertionError,
      'Query must be given');
    expect(RangerModel.hasReference({})).toBe(false);
    expect(RangerModel.hasReference({reference: ''})).toBe(false);
    expect(RangerModel.hasReference({reference: '/path/to/ref.fa'})).toBe(true);
  });
  it('provide alternate reference repository', function() {
    config.provide( () => {return {references: '/test/reference/root'};});
    let query = {reference: 'references/example/path/fasta.fa'};
    expect(RangerModel.fixReference(query)).toEqual('/test/reference/root/references/example/path/fasta.fa');

    config.provide( () => {return {};} );
  });
});

describe('Creating object instance', function() {
  beforeAll(function() {
    fse.ensureDirSync(tmpDir);
    config.provide(dummy);
  });

  afterAll(function() {
    fse.removeSync(tmpDir);
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
  let m;
  beforeAll(function() {
    fse.ensureDirSync(tmpDir);
    config.provide(dummy);
    m = new RangerModel();
  });

  afterAll(function() {
    fse.removeSync(tmpDir);
  });

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
