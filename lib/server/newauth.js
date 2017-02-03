'use strict';

const http = require('http');
const EventEmitter = require('events');
const url = require('url');

const config = require('../config.js');

const AUTH_EVENT_NAME           = 'authorised';
const AUTH_FAILED_EVENT_NAME    = 'failed';

class NewDataAccess extends EventEmitter {
  constructor() {
    super();
  }

  _authErr(e) {
    this.emit(AUTH_FAILED_EVENT_NAME,
      'Failed to get authorisation info: ' + e);
  }

  authorise(token, accessGroups) {
    let authurl = config.provide().get('authurl');
    let reqbody = {token: token, groups: accessGroups};
    let reqopts = url.parse(authurl);
    reqopts.method = 'POST';
    reqopts.path = '/checkToken';
    let req = http.request(reqopts, ( response ) => {
      let body = '';
      response.on('error', this._authErr);
      response.on('data', (chunk) => {
        body += chunk;
      });
      response.on('end', () => {
        if (response.statusCode !== 200) {
          this._authErr(new Error(
            `received ${response.statusCode} from authorisation server`));
        }
        let decision;
        try {
          decision = JSON.parse(body);
          if (decision.ok === true) {
            this.emit(AUTH_EVENT_NAME);
          } else {
            this._authErr(new Error('Not authorised for those files'));
          }
        } catch (e) {
          this._authErr(new Error('Server response was invalid JSON'));
        }
      });
    });
    req.setHeader('Content-Type', 'application/json');
    req.end(JSON.stringify(reqbody));
  }
}

module.exports = NewDataAccess;
