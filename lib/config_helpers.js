"use strict";

const assert = require('assert');
const fs     = require('fs');

let _validateAuthConfiguration = ( opts ) => {
  let authurl = opts.get('authurl');
  if ( authurl ) {
    if ( !authurl.startsWith('https') ) {
      throw Error('Only HTTPS protocol should be used for auth server');
    }

    ['auth_cert', 'auth_key'].forEach( optname => {
      assert(
        opts.get(optname),
        `'${optname}' is required when using an https Auth-server URL`
      );
    });

    ['auth_cert', 'auth_key', 'auth_ca'].forEach( optname => {
      let path = opts.get(optname);
      try {
        fs.accessSync(path, fs.R_OK);
      } catch (e) {
        throw Error(
          `File '${path}' is not readable for option '${optname}'`
        );
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

let setAuthConfiguration = ( config ) => {
  _validateAuthConfiguration(config);

  let auth_server_config = {
    url:         undefined,
    ca:          undefined,
    client_cert: undefined,
    client_key:  undefined,
    client_key_passphrase: undefined
  };

  if ( config.get('authurl') ) {
    auth_server_config.url = config.get('authurl');
    config.set('authurl', null);
  }
  if ( config.get('auth_cert') ) {
    auth_server_config.client_cert = fs.readFileSync(config.get('auth_cert'));
    config.set('auth_cert', null);
  }
  if ( config.get('auth_ca') ) {
    auth_server_config.ca = fs.readFileSync(config.get('auth_ca'));
    config.set('auth_ca', null);
  }
  if ( config.get('auth_key') ) {
    auth_server_config.client_key = fs.readFileSync(config.get('auth_key'));
    config.set('auth_key', null);
  }
  if ( config.get('auth_key_passphrase') ) {
    auth_server_config.client_key_passphrase = config.get('auth_key_passphrase');
    config.set('auth_key_passphrase', null);
  }

  auth_server_config.toMaskedString = () => {
    let tmp = JSON.parse(JSON.stringify(auth_server_config));
    tmp.client_key_passphrase = '*****';
    return JSON.stringify(tmp);
  };

  config.set('auth_server', auth_server_config);
};

module.exports = {
  setAuthConfiguration
};
