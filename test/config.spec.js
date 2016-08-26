/* globals describe, it, expect */

"use strict";
const assert  = require('assert');
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
    expect( options.get('tempdir') === '/tmp' ).toBe(true); // Provided in config.js
    expect( () => {options = config.provide( () => {return {tempdir: '/anothertmp'};} );} ).not.toThrow();
    expect( options.get('tempdir') === '/anothertmp' ).toBe(true);
  });
  it('Configs can be passed from a json file', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return { configfile: path.resolve(__dirname, 'testConfig.json') };
    });} ).not.toThrow();
    expect( options.get('testConfig') ).toBe(true);
  });
  it('Configs from function will overwrite those from json file', function() {
    let options;
    expect( () => {options = config.provide( () => {
      return {
        configfile: path.resolve(__dirname, 'testConfig.json'),
        testConfig: false
      };
    });}).not.toThrow();
    expect( options.get('testConfig') ).toBe(false);
  });
});
