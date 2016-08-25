/* globals describe, it, expect, afterEach */

"use strict";
const assert  = require('assert');
const path    = require('path');
const config  = require('../lib/config.js');

describe('Building options', function() {
  afterEach(function() {
   config.flush();
  }); 
  it('Function must be passed', function() {
    expect( () => {config.build();} ).toThrowError(assert.AssertionError,
      'generateConfigs must be a function');
  });
  it('Function must return an object', function() {
    expect( () => {config.build( () => {return 1;} );} ).toThrowError(assert.AssertionError,
      'generateConfigs must return an object');
    expect( () => {config.build( () => {return {};} );} ).not.toThrow();
  });
  it('Can get options passed via function', function() {
    let options;
    expect( () => {options = config.build( () => {return {test: 'pass'};} );} ).not.toThrow();
    expect( options.get('test') === 'pass' ).toBe(true);
  });
  it('Options passed via function will overwrite defaults', function() {
    let options;
    expect( () => {options = config.build( () => {return {};} );} ).not.toThrow();
    expect( options.get('tempdir') === '/tmp' ).toBe(true); // Provided in config.js
    config.flush();
    expect( () => {options = config.build( () => {return {tempdir: '/anothertmp'};} );} ).not.toThrow();
    expect( options.get('tempdir') === '/anothertmp' ).toBe(true);
  });
  it('Configs can be passed from a json file', function() {
    let options;
    expect( () => {options = config.build( () => {
      return { configfile: path.resolve(__dirname, 'testConfig.json') };
    });} ).not.toThrow();
    expect( options.get('testConfig') ).toBe(true);
  });
  it('Configs from function will overwrite those from json file', function() {
    let options;
    expect( () => {options = config.build( () => {
      return {
        configfile: path.resolve(__dirname, 'testConfig.json'),
        testConfig: false
      };
    });}).not.toThrow();
    expect( options.get('testConfig') ).toBe(false);
  });
});
