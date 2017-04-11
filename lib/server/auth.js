'use strict';

/**
 * auth module.
 * @module server/auth
 *
 * @description Authorisation procedure.
 *              Makes requests to the /checkToken and /checkUser endpoints
 *              on a npg_sentry server (found via the authurl configuration
 *              option) to determine whether a user (identified by either
 *              their email address when authtype === 'authuser' or a token
 *              when authtype === 'authtoken') is a member of every group
 *              passed.
 *
 * @example <caption>Example usage of the auth module.</caption>
 *   const DataAccess = require('../lib/server/auth.js');
 *
 *   // Create a new object.
 *   var authtype = 'authuser' || 'authtoken';
 *   var da = new DataAccess(authtype);
 *   // Register listeners for events this object emits.
 *   da.on('authorised', (username) => {
 *          console.log(`User ${username} is given access`);
 *        });
 *   da.on('failed', (username, message) => {
 *          console.log(
 *            `Authorisation failed for user '${username}': ${message}`);
 *        });
 *   // Authorise a user for access for data with data access group ids
     // given by the second attribute. Group ids should be strings.
 *   da.authorise('alice', ["6", "34"]);
 *
 * @author Marina Gourtovaia
 * @author Andrew Nowak
 * @copyright Genome Research Limited 2016
 */

const EventEmitter = require('events');
const url    = require('url');

const LOGGER  = require('winston');
const request = require('request');

const config    = require('../config.js');
const constants = require('../constants.js');

const AUTH_EVENT_NAME        = 'authorised';
const AUTH_FAILED_EVENT_NAME = 'failed';

const DEFAULT_AUTH_ERROR_MSG = 'Unexpected error while processing authorisation';

/** Class encupsulating authorisation. */
class DataAccess extends EventEmitter {

  /**
   * Creates an DataAccess type object.
   * @param authtype - string signifying how to authorise user
   */
  constructor(authtype) {
    super();
    if (!authtype) {
      throw new ReferenceError('Authorisation type is required');
    }
    if ([constants.AUTH_TYPE_USER, constants.AUTH_TYPE_TOKEN].indexOf(authtype) < 0) {
      throw new ReferenceError('Unknown authorisation type');
    }
    this.authtype = authtype;
  }

  _authErr(error, errorMessage) {
    let message = errorMessage || DEFAULT_AUTH_ERROR_MSG;
    if ( error !== null ) {
      LOGGER.error('Failed to get authorisation info: ' + error);
    }
    this.emit(AUTH_FAILED_EVENT_NAME, message);
  }

  _getName(uname) {
    uname = uname.trim();
    if (uname) {
      let domain = config.provide().get('emaildomain');
      if (domain) {
        // Allow alphanumerics, _, - and . in a username.
        // Ensure special characters are properly escaped.
        let re = new RegExp(
          '([\\w\\.-]+)@' +
          domain.replace(/\./g, '\\.')
          + '$', 'im');
        let result = uname.match(re);
        uname = result ? result[1] : null;
      }
    }
    return uname;
  }

  _validateObjHasOnly(obj, expectedProps) {
    for (let i = 0; i < expectedProps.length; i++) {
      let property = expectedProps[i];
      if ( !(property in obj) ) {
        this._authErr(new Error(
          `Auth server response must contain value for property: '${property}'`,
          DEFAULT_AUTH_ERROR_MSG
        ));
        return false;
      }
    }
    for ( let property in obj ) {
      if ( obj.hasOwnProperty(property) && expectedProps.indexOf(property) === -1 ) {
        this._authErr(new Error(
          `Auth server response contains unexpected value: '${property}'`,
          DEFAULT_AUTH_ERROR_MSG
        ));
        return false;
      }
    }
    return true;
  }

  /**
   * Authorise the user and emit an appropriate event.
   * @param identifier - username or token
   * @param accessGroups - an array of access groups to authorise
   *                       the user against
   */
  authorise(identifier, accessGroups) {
    if (!identifier) {
      throw new ReferenceError('Identifier is required');
    }
    if (!accessGroups || (accessGroups instanceof Array === false) ||
       (accessGroups.length === 0)) {
      throw new Error('Access groups array is not available');
    }

    if (!accessGroups.every( (id) => { return id; })) {
      this._authErr(new Error('Some access group ids are not defined'));
      return;
    }

    let agroupIds = accessGroups.slice();
    agroupIds.sort();
    agroupIds = agroupIds.filter( function(item, index, thisArray) {
      return (index === 0) ? 1 : ((item === thisArray[index - 1]) ? 0 : 1);
    });
    LOGGER.debug('ACCESS GROUP IDS: ' + agroupIds.join(' '));

    let auth_server = config.provide().get('auth_server');
    let authurl = auth_server.url;
    let reqbody = {groups: agroupIds};
    let reqopts = {};
    reqopts.url     = url.parse(authurl);
    reqopts.method  = 'POST';
    reqopts.headers = {
      'Content-Type': 'application/json'
    };

    if ( authurl.startsWith(constants.SECURE_PROTOCOL) ) {
      let auth_cert = auth_server.client_cert;
      if ( auth_cert ) {
        reqopts.cert = auth_cert;
      }
      let auth_key = auth_server.client_key;
      if ( auth_key ) {
        reqopts.key = auth_key;
        let auth_passphrase = auth_server.client_key_passphrase;
        if ( auth_passphrase ) {
          reqopts.passphrase = auth_passphrase;
        }
      }
      let auth_ca = auth_server.ca;
      if ( auth_ca ) {
        reqopts.ca = auth_ca;
      }
    } else { // This should never happen because of validation in config
      LOGGER.error('Auth server can only run as https service');
      process.exit(1);
    }

    if (this.authtype === constants.AUTH_TYPE_TOKEN) {
      reqbody.token = identifier.trim();
      reqopts.url.path = constants.AUTH_URL_TOKEN;
    } else if (this.authtype === constants.AUTH_TYPE_USER) {
      reqbody.user = this._getName(identifier);
      reqopts.url.path = constants.AUTH_URL_USER;
    }
    if (!(reqbody.token || reqbody.user)) {
      return this._authErr(new Error(`Invalid identifier "${identifier}"`));
    }

    reqopts.body = JSON.stringify(reqbody);

    request(reqopts, ( error, response, body ) => {
      if ( error ) {
        this._authErr( error );
      } else {
        if ( response.statusCode !== 200 ) {
          this._authErr(
            new Error(`received ${response.statusCode} from authorisation server`),
            DEFAULT_AUTH_ERROR_MSG
          );
        } else {
          let decision;
          try {
            decision = JSON.parse(body);
          } catch (e) {
            this._authErr(
              new Error(`Server response was invalid JSON ${e}`),
              DEFAULT_AUTH_ERROR_MSG
            );
            return;
          }

          // Check only expected values in body
          let expectedProperties = ['ok'];
          if ( !this._validateObjHasOnly(decision, expectedProperties) ) {
            return;
          }

          if ( typeof decision.ok !== 'boolean' ) {
            return this._authErr(new Error(
              'Auth server response does not match expected format',
              DEFAULT_AUTH_ERROR_MSG
            ));
          }

          if (decision.ok === true) {
            this.emit(AUTH_EVENT_NAME);
          } else {
            this._authErr(
              null,
              'Not authorised for those files'
            );
          }
        }
      }
    });
  }
}

module.exports = DataAccess;
