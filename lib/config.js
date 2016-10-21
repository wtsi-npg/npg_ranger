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
const urlParser   = require('url');
const assert      = require('assert');
const moment      = require('moment');
const GetOpt      = require('node-getopt');

const DEFAULT_PROTOCOL = 'http:';
const RO_OPTION_KEY    = 'config_ro';

const optionsList = [
  ['p','port=PORT'        ,'PORT or socket which server listens on'],
  ['m','mongourl=URI'     ,'URI to contact mongodb'],
  ['t','tempdir=PATH'     ,'PATH of temporary directory'],
  ['H','hostname=HOST'    ,'override hostname with HOST'],
  ['' ,'multiref'         ,'allow merging files with differing references. disables VCF output'],
  ['c','configfile=PATH'  ,'PATH of configfile'],
  ['g','timeout=SECONDS'  ,'SECONDS to wait before killing child processes'],
  ['r','references=PATH'  ,'PATH to dir containing reference repository'],
  ['s','skipauth'         ,'skip authorisation steps'],
  ['' ,'anyorigin'        ,'accept CORS requests from any origin'],
  ['d','debug'            ,'debugging mode for this server'],
  ['V','version'          ,'print version and exit'],
  ['h','help'             ,'display this help']
];

const defaultOptions = {
  config_ro:  false,
  originlist: null,
  proxylist:  null,
  protocol:   DEFAULT_PROTOCOL,
  hostname:   os.hostname() || 'localhost',
  mongourl:   'mongodb://localhost:27017/imetacache',
  timeout:    3,
  mongoopt: {
    db: {
      numberofRetries: 5
    },
    server: {
      auto_reconnect: true,
      poolSize:       40,
      socketOptions:  {
        connectTimeoutMS: 5000
      }
    },
    replSet: {},
    mongos: {}
  }
};

var options;

var fromCommandLine = () => {
  var opt = new GetOpt(optionsList).bindHelp().parseSystem();
  return opt.options;
};

/**
 * <p>Builds and/or retrieves the options object. Validates options after they are built
 * If an error happens during validation, the options preceeding this reset attempt
 * are retained.</p>
 * <p><strong>WARNING</strong>: if a function is provided, all configurations
 * originating from a previous function <strong>WILL BE LOST</strong>, including
 * those from a config file (unless the same config file is provided in the new
 * function).</p>
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
 * // Make the configuration read-only
 * options = config.provide( () = {
 *   return {baz: true};
 * }, true); // Makes the configuration read-only
 * options = config.provide( () = {
 *   return {baz: false};
 * }); // Will trow Error
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
 * @param {Boolean} immutable - Determines if the resulting options object will
 *                              be immutable.
 *
 * @throws {AssertionError} Provided parameter was not a function that returned
 *                          an object, or no parameter was provided but options
 *                          had not yet been initialized
 * @throws {RangerError}    If an error in validating options happens
 * @throws {Error}          If an attempt to overwrite a read-only configuration
 *                          is made
 *
 * @return {object} Object which can be queried with <code>.get('key')</code> to
 *                  find the values of a given setting.
 */
var provide = function(generateConfigs, immutable) {
  if (generateConfigs) {
    if (options && options.get(RO_OPTION_KEY)) {
      throw new Error('Attempt to overwrite original configuration');
    }
    assert(typeof generateConfigs === 'function', 'parameter must be a function');
    let opts = generateConfigs();
    assert(typeof opts === 'object', 'parameter must return an object');
    let tmpOptions = new ConfigChain(
      opts,
      opts.configfile
        ? path.resolve(opts.configfile)
        : null,
      defaultOptions
    );
    if (typeof immutable !== 'undefined') {
      assert(typeof immutable === 'boolean', 'immutable must be boolean');
      tmpOptions.set('config_ro', immutable);
    }
    adjustOptions(tmpOptions);
    if (tmpOptions.get(RO_OPTION_KEY)) {
      tmpOptions.set = () => {
        throw new Error('Attempt to change read-only configuration');
      };
    }
    options = tmpOptions;
  }
  assert(options, 'Options is undefined');
  return options;
};

function adjustOptions(opts) {

  let validateURL = (u, name) => {
    u = u.trim();
    if (!u) {
      throw new RangeError(`Empty string in '${name}'`);
    }
    try {
      let urlObj = urlParser.parse(u);
      if (!urlObj.protocol) {
        throw new Error('Protocol is absent');
      }
      if (urlObj.protocol !== opts.get('protocol')) {
        throw new Error('URL protocol should match server protocol');
      }
      if (!(urlObj.host || urlObj.hostname)) {
        throw new Error('Server host is absent');
      }
      if (urlObj.pathname && urlObj.path && urlObj.pathname.length > 1) {
        throw new Error('Path cannot be present');
      }
      if (urlObj.search) {
        throw new Error('Search string cannot be present');
      }
      if (urlObj.hash) {
        throw new Error('Hash tag cannot be present');
      }
      u = urlParser.format(urlObj);
      u = u.replace(/\/*$/, ''); // Drop trailing slash
    } catch (e) {
      throw new RangeError(`Invalid URL ${u} in ${name}: ${e}`);
    }
    return u;
  };

  let setTempDir = () => {
    opts.set('tempdir', (opts.get('tempdir') || tempFilePath('npg_ranger_')));
  };

  let setPort = () => {
    if (!opts.get('port')) {
      opts.set('port', path.join(opts.get('tempdir'), 'npg_ranger.sock'));
    }
    if (Number.isInteger(Number.parseInt(opts.get('port')))) {
      opts.set('port', Number.parseInt(opts.get('port')));
    }
  };

  let setACAOrigin = () => {
    let or = opts.get('originlist');
    if (opts.get('anyorigin')) {
      if (typeof opts.get('anyorigin') != 'boolean') {
        throw new RangeError("'anyorigin' should be a boolean type value");
      }
      if (or) {
        throw new RangeError("'anyorigin' and 'originlist' options cannot both be set");
      }
      if (!opts.get('skipauth')) {
        throw new RangeError("'anyorigin' option cannot be set if authorization is performed");
      }
    } else if (or) {
      assert(or instanceof Array, "'originlist' should be an array");
      opts.set('originlist',
        or.length === 0
        ? null
        : or.filter((el) => {return el;})
            .map((el) => {
              return validateURL(el, 'originlist');
            })
      );
    }
  };

  let setProxies = () => {
    let proxies = opts.get('proxylist');
    if (proxies) {
      assert(proxies instanceof Object, "'proxilist' should be a hash");
      let urls = Object.keys(proxies);
      if (urls.length === 0) {
        proxies = null;
      } else {
        urls.forEach((el) => {
          if (el) {
            let u = validateURL(el, 'proxylist');
            if (u !== el) {
              proxies[u] = proxies[el];
              delete proxies[el];
            }
          } else {
            throw new RangeError('Empty or zero url in proxilist');
          }
        });
      }
      opts.set('proxylist', proxies);
    }
  };

  setTempDir();
  setPort();
  setACAOrigin();
  setProxies();
}

var _formatDate = () => {
  return moment().format('YYYYMMDD_HHmmssSS');
};

var tempFilePath = function(prefix) {
  prefix = prefix ? prefix : '';
  return path.join(os.tmpdir(), prefix + _formatDate());
};

var logOpts = function() {
  let loIndex = 1;

  let names = optionsList
    .filter(function(el) {return !el[loIndex].startsWith('help');})
    .sort(function(el1, el2) {return el1[loIndex].localeCompare(el2[loIndex]);})
    .map(function(el) {
      let desc = el[loIndex];
      return desc.split('=', 1).join('');
    });
  // Push all defaultOptions into names
  Array.prototype.push.apply(names, Object.keys(defaultOptions).sort());

  return "\n" + names
    .map(function(name) {
      return name + '=' + JSON.stringify(options.get(name));
    }).join("\n");
};

module.exports = {
  fromCommandLine: fromCommandLine,
  provide:         provide,
  tempFilePath:    tempFilePath,
  logOpts:         logOpts
};
