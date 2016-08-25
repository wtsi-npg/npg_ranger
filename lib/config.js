"use strict";

/**
 * configuration module.
 * @module config.js
 *
 * @description Configuration module for retrieving and providing settings.
 *
 * @example <caption>Example usage.</caption>
 *   const config = require('../lib/config.js');
 *   // Build the config store.
 *   // func should be a function returning an object.
 *   // The values from this object will override any
 *   // found elsewhere, in config files or similar.
 *   config.build(func);
 *   // Retrieve the 'debug' setting
 *   options.get('debug');
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

  return opt;
};

var options;

var build = function( generateConfigs ) {
  if ( !options ) {
    assert(typeof generateConfigs === 'function', 'generateConfigs must be a function');
    let opts = generateConfigs();
    assert(typeof opts === 'object', 'generateConfigs must return an object');
    console.log('defining options');
    options = new ConfigChain(
      opts,
      opts.configfile
        ? path.resolve(__dirname, opts.configfile)
        : null,
      { // TODO maybe remove
        port: path.join(os.tmpdir(), process.env.USER, 'npg_ranger_sock'),
        hostname: os.hostname(),
        tempdir: '/tmp',
        mongourl: 'mongodb://sf2-farm-srv1:27017/imetacache',
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
  return options;
};

var flush = function() {
  options = undefined;
};

module.exports = {
  fromCommandLine: fromCommandLine,
  build: build,
  flush: flush
};

/*
const PORT               = opt.options.port || opt.argv[0]
                           || path.join(os.tmpdir(), process.env.USER, 'npg_ranger.sock');
const HOST               = opt.options.hostname || os.hostname() || 'localhost';
const MONGO              = opt.options.mongourl || 'mongodb://sf2-farm-srv1:27017/imetacache';
const TEMP_DATA_DIR_NAME = 'npg_ranger_data';
const TEMP_DATA_DIR      = opt.options.tempdir || path.join(os.tmpdir(), process.env.USER, TEMP_DATA_DIR_NAME);
*/
