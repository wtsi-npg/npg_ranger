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
 *   config.provide(func);
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
const moment      = require('moment');
const GetOpt      = require('node-getopt');

var optionsList = [
      ['p','port=PORT'        ,'PORT or socket which server listens on'],
      ['m','mongourl=URI'     ,'URI to contact mongodb'],
      ['t','tempdir=PATH'     ,'PATH of temporary directory'],
      ['H','hostname=HOST'    ,'override hostname with HOST'],
      ['' ,'multiref'         ,'allow merging files with differing references. disables VCF output'],
      ['c','configfile=PATH'  ,'PATH of configfile'],
      ['g','timeout=SECONDS'  ,'SECONDS to wait before killing child processes'],
      ['r','references=PATH'  ,'PATH to dir containing reference repository'],
      ['s','skipauth'         ,'skip authorisation steps'],
      ['' ,'cors'             ,
        'add CORS headers (will work only if both skipauth is true and not using a socket)'],
      ['d','debug'            ,'debugging mode for this server'],
      ['h','help'             ,'display this help']
];

var options;

var fromCommandLine = () => {
  var opt = new GetOpt(optionsList).bindHelp().parseSystem();
  return opt.options;
};

/**
 * <p>Builds and/or retrieves the options object.</p>
 * <p>WARNING: if a function is provided, all configurations originating from a previous
 * function WILL BE LOST, including those from a config file (unless the same config file
 * is provided in the new function).</p>
 *
 * @example
 * var options;
 * // Build and retrieve the options object.
 * function foo() {
 *   return {baz: true};
 * }
 * options = config.provide(foo);
 * options.get('baz'); // true
 * // Retrieve the options object without clearing the options object.
 * options = config.provide();
 * options.get('baz'); // true
 * // Clear the options object and rebuild it with the new function
 * function bar() {
 *   return {baz: false};
 * }
 * options = config.provide(bar);
 * options.get('baz'); // false
 *
 *
 * @param {function} generateConfigs - Function that returns an object. The keys
 *                                     and values of this object will be used to
 *                                     as the keys and values of the resulting
 *                                     options object.
 *                                     <br>
 *                                     If this function is not supplied,
 *                                     <code>provide</code> will attempt to
 *                                     retrieve an already created options
 *                                     object.
 *                                     <br>
 *                                     Configs can be retrieved from a json file if the
 *                                     returned object provides a path specified by the
 *                                     'configfile' key, but the returned object will
 *                                     always take precedence.
 *
 * @throws {AssertionError} Provided parameter was not a function that returned
 *                          an object, or no parameter was provided but options
 *                          had not yet been initialized.
 *
 * @return {object} Object which can be queried with <code>.get('key')</code> to
 *                  find the values of a given setting.
 */
var provide = function( generateConfigs ) {
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
      {
        hostname: os.hostname() || 'localhost',
        mongourl: 'mongodb://localhost:27017/imetacache',
        timeout: 3,
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
    options.set('tempdir', (options.get('tempdir') || tempFilePath('npg_ranger_')));
    if (!options.get('port')) {
      options.set('port', path.join(options.get('tempdir'), 'npg_ranger.sock'));
    }
  }
  assert(options, 'Options is undefined');
  return options;
};

var _formatDate = () => {
  return moment().format('YYYYMMDD_HHmmssSS');
};

var tempFilePath = function(prefix) {
  prefix = prefix ? prefix : '';
  return path.join(os.tmpdir(), prefix + _formatDate());
};

var logOpts = function() {
  let loIndex = 1;
  return optionsList
    .filter(function(el) {return !el[loIndex].startsWith('help');})
    .sort(function(el1, el2) {return el1[loIndex].localeCompare(el2[loIndex]);})
    .map(function(el) {
      let desc = el[loIndex];
      let longDescription = desc.split('=', 1).join('');
      return longDescription + '=' + options.get(longDescription);
    }).join(', ');
};

module.exports = {
  fromCommandLine: fromCommandLine,
  provide:         provide,
  tempFilePath:    tempFilePath,
  logOpts:         logOpts
};
