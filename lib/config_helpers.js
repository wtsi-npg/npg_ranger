"use strict";

const assert = require('assert');
const fs     = require('fs');

const constants = require('./constants');

const MASK = '*****';

let _validateAuthConfiguration = ( opts ) => {
  let authurl = opts.get('authurl');
  if ( authurl ) {
    if ( !authurl.startsWith(constants.SECURE_PROTOCOL) ) {
      throw Error('Only HTTPS protocol should be used for auth server');
    }

    if ( opts.get('auth_cert') || opts.get('auth_key') ) {
      ['auth_cert', 'auth_key'].forEach( optname => {
        assert(
          opts.get(optname),
          `'${optname}' is required when using auth-cert/key pair`
        );
      });
    }

    ['auth_cert', 'auth_key', 'auth_ca'].forEach( optname => {
      let path = opts.get(optname);
      if ( path ) {
        try {
          fs.accessSync(path, fs.R_OK);
        } catch (e) {
          throw Error(
            `File '${path}' is not readable for option '${optname}'`
          );
        }
      }
    });
  } else {
    ['auth_cert',
     'auth_key',
     'auth_ca',
     'auth_key_passphrase'].forEach( (name) => {
      if (opts.get(name)) {
        throw Error(
          `A value for option '${name}' is available but no authurl was provided`
        );
      }
    });
  }
};

/**
 * Wrap all auth server configuration into a single object and store it in the
 * original configuration under <code>auth-server</code> key.
 * @param {object} opts <i>config-chain</i> configuration object
 */
let setAuthConfiguration = ( opts ) => {
  _validateAuthConfiguration(opts);

  let auth_server_config = {
    url:         undefined,
    ca:          undefined,
    client_cert: undefined,
    client_key:  undefined,
    client_key_passphrase: undefined
  };

  if ( opts.get('authurl') ) {
    auth_server_config.url = opts.get('authurl');
    opts.set('authurl', null);
  }
  if ( opts.get('auth_cert') ) {
    auth_server_config.client_cert = fs.readFileSync(opts.get('auth_cert'));
    opts.set('auth_cert', null);
  }
  if ( opts.get('auth_ca') ) {
    auth_server_config.ca = fs.readFileSync(opts.get('auth_ca'));
    opts.set('auth_ca', null);
  }
  if ( opts.get('auth_key') ) {
    auth_server_config.client_key = fs.readFileSync(opts.get('auth_key'));
    opts.set('auth_key', null);
  }
  if ( opts.get('auth_key_passphrase') ) {
    auth_server_config.client_key_passphrase = opts.get('auth_key_passphrase');
    opts.set('auth_key_passphrase', null);
  }

  auth_server_config.toMaskedString = () => {
    let tmp = JSON.parse(JSON.stringify(auth_server_config)); // clone
    if ( tmp.client_key_passphrase ) {
      tmp.client_key_passphrase = MASK;
    }
    return JSON.stringify(tmp);
  };

  opts.set('auth_server', auth_server_config);
};

/**
 * Wrap all ssl server configuration into a single object and store it in the
 * original configuration under <code>server_ssl</code> key.
 * @param {object} opts <i>config-chain</i> configuration object
 */
let setSSLConfiguration = (opts) => {
  let startssl = opts.get('startssl');
  if ( startssl && opts.get('protocol') != constants.SECURE_PROTOCOL ) {
    opts.set('protocol', constants.SECURE_PROTOCOL);
  }

  ['secure_key', 'secure_passphrase', 'secure_cert'].forEach( (optname) => {
    if ( startssl ) {
      if ( optname !== 'secure_passphrase' ) { // passphrase is optional
        let path = opts.get(optname);
        assert(path, `'${optname}' is required when using 'startssl' option`);
        try {
          fs.accessSync(path, fs.R_OK);
        } catch (e) {
          throw Error(
            `File '${path}' is not readable for option '${optname}'`
          );
        }
      }
    } else {
      if ( opts.get(optname) ) {
        throw new RangeError(`'${optname}' option requires startssl to be true`);
      }
    }
  });

  let ssl_config = {
    server_cert: undefined,
    server_key:  undefined,
    server_key_passphrase: undefined
  };

  if ( opts.get('secure_key') ) {
    ssl_config.server_key = opts.get('secure_key');
    opts.set('secure_key', null);
  }
  if ( opts.get('secure_cert') ) {
    ssl_config.server_cert = opts.get('secure_cert');
    opts.set('secure_cert', null);
  }
  if ( opts.get('secure_passphrase') ) {
    ssl_config.server_key_passphrase = opts.get('secure_passphrase');
    opts.set('secure_passphrase', null);
  }

  ssl_config.toMaskedString = () => {
    let tmp = JSON.parse(JSON.stringify(ssl_config)); // clone
    if ( tmp.server_key_passphrase ) {
      tmp.server_key_passphrase = MASK;
    }
    return JSON.stringify(tmp);
  };

  opts.set('server_ssl', ssl_config);
};

module.exports = {
  setAuthConfiguration,
  setSSLConfiguration
};
