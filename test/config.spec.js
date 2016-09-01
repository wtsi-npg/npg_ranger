/* globals describe, it, expect */

"use strict";
const assert  = require('assert');
const os      = require('os');
const path    = require('path');
const config  = require('../lib/config.js');

describe('Building options', function() {
  it('Must define options before retrieval', function() {
    expect( () => {config.provide();} ).toThrowError(assert.AssertionError,
      'Options is undefined');
  });
  it('If there is a parameter, it must be a function', function() {
    expect( () => {config.provide( 1 );} ).toThrowError(assert.AssertionError,
      'parameter must be a function');
  });
  it('Function must return an object', function() {
    expect( () => {config.provide( () => {return 1;} );} ).toThrowError(assert.AssertionError,
      'parameter must return an object');
    expect( () => {config.provide( () => {return {};} );} ).not.toThrow();
  });
  it('Can get options passed via function', function() {
    let options;
    expect( () => {options = config.provide( () => {return {test: 'pass'};} );} ).not.toThrow();
    expect( options.get('test') === 'pass' ).toBe(true);
  });
  it('Options passed via function will overwrite defaults', function() {
    let options;
    expect( () => {options = config.provide( () => {return {};} );} ).not.toThrow();
    expect( options.get('mongourl') === 'mongodb://localhost:27017/imetacache' ).toBe(true); // Provided in config.js
    expect( () => {options = config.provide( () => {return {mongourl: 'newmongourl'};} );} ).not.toThrow();
    expect( options.get('mongourl') === 'newmongourl' ).toBe(true);
  });
  it('Configs can be passed from a json file', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return { configfile: path.resolve(__dirname, 'server', 'data', 'testConfig.json') };
    });} ).not.toThrow();
    expect( options.get('testConfig') ).toBe(true);
  });
  it('Configs from function will overwrite those from json file', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return {
        configfile: path.resolve(__dirname, 'server', 'data', 'testConfig.json'),
        testConfig: false
      };
    });}).not.toThrow();
    expect( options.get('testConfig') ).toBe(false);
  });
  it('tempdir, port unspecified, defaults created', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return {};
    });}).not.toThrow();
    expect( options.get('tempdir').startsWith(path.join(os.tmpdir(), 'npg_ranger_')) ).toBe(true);
    expect( options.get('port').startsWith(path.join(os.tmpdir(), 'npg_ranger_')) ).toBe(true);
    expect( options.get('port').endsWith('npg_ranger.sock') ).toBe(true);
  });
  it('tempdir specified, port default', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return {tempdir: config.tempFilePath('npg_ranger_config_test_')};
    });}).not.toThrow();
    expect( options.get('tempdir').startsWith(path.join(os.tmpdir(), 'npg_ranger_config_test')) ).toBe(true);
    expect( options.get('port').startsWith(path.join(os.tmpdir(), 'npg_ranger_config_test'))).toBe(true);
    expect( options.get('port').endsWith('npg_ranger.sock') ).toBe(true);
  });
  it('tempdir default, port specified', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return {port: '45678'};
    });}).not.toThrow();
    expect( options.get('tempdir').startsWith(path.join(os.tmpdir(), 'npg_ranger_')) ).toBe(true);
    expect( options.get('port') == '45678' ).toBe(true);
  });
  it('tempdir specified, port specified', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return {
        tempdir: config.tempFilePath('npg_ranger_config_test_'),
        port:    '45678'
      };
    });}).not.toThrow();
    expect( options.get('tempdir').startsWith(path.join(os.tmpdir(), 'npg_ranger_config_test_')) ).toBe(true);
    expect( options.get('port') == '45678' ).toBe(true);
  });
});

describe('Creating temp file path', function() {
  it('Without prefix', function() {
    let temppath = config.tempFilePath();
    expect( temppath.startsWith(os.tmpdir()) ).toBe(true);
    expect( temppath ).toMatch(/\/\d*$/);
  });
  it('With prefix', function() {
    let temppath = config.tempFilePath('npg_ranger_config_test_');
    expect( temppath.startsWith(path.join(os.tmpdir(), 'npg_ranger_config_test_')) ).toBe(true);
    expect( temppath ).toMatch(/\/npg_ranger_config_test_\d*$/);
  });
});
