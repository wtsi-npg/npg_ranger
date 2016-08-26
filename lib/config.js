"use strict";

/**
 * configuration module.
 * @module config
 *
 * @description Configuration module for retrieving and providing settings.
 *
 * @example <caption>Example usage.</caption>
 *   const config = require('../lib/config.js');
 *   // Build the config store.
 *   // func should be a function returning an object.
 *   // The values from this object will override any
 *   // found elsewhere, in config files or similar.
 *   function func() {
 *     return { debug: true };
 *   }
 *   config.build(func);
 *   // Retrieve the 'debug' setting
 *   options.get('debug'); // true
 *
 * @author Andrew Nowak
 * @copyright Genome Research Limited 2016
 */

const ConfigChain = require('config-chain');
const path        = require('path');
const os          = require('os');
const assert      = require('assert');

const GetOpt  = require('node-getopt');

var fromCommandLine = () => {
  var opt = new GetOpt([
      ['p','port=PORT'        ,'PORT or socket which server listens on'],
      ['m','mongourl=URI'     ,'URI to contact mongodb'],
      ['t','tempdir=PATH'     ,'PATH of temporary directory'],
      ['H','hostname=HOST'    ,'override hostname with HOST'],
      ['c','configfile=PATH'  ,'PATH of configfile'],
      ['s','skipauth'         ,'skip authorisation steps'],
      ['d','debug'            ,'debugging mode for this server'],
      ['h','help'             ,'display this help']
  ]).bindHelp().parseSystem();

  return opt.options;
};

var options;

/**
 * <p>Builds and/or retrieves the options object.</p>
 *
 * @example
 * var options;
 * // Build and retrieve the options object.
 * function foo() {
 *   return {baz: true};
 * }
 * options = config.build(foo);
 * options.get('baz'); // true
 * // Retrieve the options object without clearing the options object.
 * options = config.build();
 * options.get('baz'); // true
 * // Clear the options object and rebuild it with the new function
 * function bar() {
 *   return {baz: false};
 * }
 * options = config.build(bar);
 * options.get('baz'); // false
 *
 *
 * @param {function} generateConfigs - Function that returns an object. The keys
 *                                     and values of this object will be used to
 *                                     as the keys and values of the resulting
 *                                     options object.
 *                                     <br>
 *                                     If this function is not supplied,
 *                                     <code>build</code> will attempt to
 *                                     retrieve an already created options
 *                                     object. If none exist, <code>build</code>
 *                                     will throw an AssertionError.
 *
 * @throws {AssertionError} Provided parameter was not a function that returned
 *                          an object, or no parameter was provided but options
 *                          had not yet been initialized.
 *
 * @return {object} Object which can be queried with <code>.get('key')</code> to
 *                  find the values of a given setting.
 */
var build = function( generateConfigs ) {
  if ( generateConfigs ) {
    assert(typeof generateConfigs === 'function', 'parameter must be a function');
    let opts = generateConfigs();
    assert(typeof opts === 'object', 'parameter must return an object');
    console.log('defining options');
    options = new ConfigChain(
      opts,
      opts.configfile
        ? path.resolve(opts.configfile)
        : null,
      { // TODO maybe remove
        port: path.join(os.tmpdir(), process.env.USER, 'npg_ranger_sock'),
        hostname: os.hostname() || 'localhost',
        tempdir: '/tmp',
        mongourl: 'mongodb://localhost:27017/imetacache',
        mongoopt: {
          db: {
            numberofRetries: 5
          },
          server: {
            auto_reconnect: true,
            poolSize: 40,
            socketOptions: {
              connectTimeoutMS: 5000
            }
          },
          replSet: {},
          mongos: {}
        }
      }
    );
  }
  assert(options, 'Options is undefined');
  return options;
};
// TODO possibly rename build? because it is not always building now

module.exports = {
  fromCommandLine: fromCommandLine,
  build: build
};
