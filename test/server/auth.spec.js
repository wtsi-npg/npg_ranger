/* globals describe, it, expect, beforeAll, afterAll, fail */

'use strict';

const http        = require('http');
const url         = require('url');

const DataAccess  = require('../../lib/server/auth.js');
const config      = require('../../lib/config.js');
const constants   = require('../../lib/constants.js');

describe('Authorisation', function() {
  let closeServ;
  let dummyAuthServ;
  var noemail_conf = function() {
    return {
      emaildomain: null,
      authurl: 'http://localhost:9898/'
    };
  };
  var email_conf = function() {
    return {
      emaildomain: 'boom.co.uk',
      authurl: 'http://localhost:9898/'
    };
  };

  var isSubset = function(arr1, arr2) {
    return arr1.every(function(val) {
      return arr2.indexOf(val) > -1;
    });
  };

  beforeAll( () => {
    config.provide(noemail_conf);

    // Provide a minimal standin for npg_sentry.
    // Authorises user 'alice' or token 'abc' for groups 1, 2 and 3.
    dummyAuthServ = http.createServer((req, res) => {
      if (req.method !== 'POST') {
        res.statusCode = 404;
        return res.end();
      }
      let body = '';
      req.on('data', (chunk) => {
        body += chunk;
      });
      let requrl = url.parse(req.url);
      let ok;
      req.on('end', () => {
        let reqdata = JSON.parse(body);
        if (requrl.pathname === constants.AUTH_URL_TOKEN) {
          if (reqdata.token === 'abc' && isSubset(reqdata.groups, ['1', '2', '3'])) {
            ok = true;
          } else {
            ok = false;
          }
          res.end(JSON.stringify({ok: ok}));
        } else if (requrl.pathname === constants.AUTH_URL_USER) {
          if (reqdata.user === 'alice' && isSubset(reqdata.groups, ['1', '2', '3'])) {
            ok = true;
          } else {
            ok = false;
          }
          res.end(JSON.stringify({ok: ok}));
        } else {
          res.statusCode = 404;
          res.end();
        }
      });
    });
    dummyAuthServ.listen(9898);
    closeServ = dummyAuthServ.close;
  });

  it('Input validation', function() {
    expect( () => {new DataAccess();} ).toThrowError(ReferenceError,
      'Authorisation type is required');
    expect( () => {new DataAccess('badauth');} ).toThrowError(ReferenceError,
      'Unknown authorisation type');
    expect( () => {new DataAccess(constants.AUTH_TYPE_USER);} ).not.toThrow();
    expect( () => {new DataAccess(constants.AUTH_TYPE_TOKEN);} ).not.toThrow();
    expect( () => {
      let da = new DataAccess(constants.AUTH_TYPE_USER);
      da.authorise();
    }).toThrowError(ReferenceError, 'Identifier is required');
    expect( () => {
      let da = new DataAccess(constants.AUTH_TYPE_TOKEN);
      da.authorise();
    }).toThrowError(ReferenceError, 'Identifier is required');
    expect( () => {
      let da = new DataAccess(constants.AUTH_TYPE_USER);
      da.authorise('alice');
    }).toThrowError(Error, 'Access groups array is not available');
    expect( () => {
      let da = new DataAccess(constants.AUTH_TYPE_USER);
      da.authorise('alice', 2);
    }).toThrowError(Error, 'Access groups array is not available');
    expect( () => {
      let da = new DataAccess(constants.AUTH_TYPE_USER);
      da.authorise('alice', []);
    }).toThrowError(Error, 'Access groups array is not available');
  });

  it('Authorisation failed - username is all whitespace', function(done) {
    var da = new DataAccess(constants.AUTH_TYPE_USER);
    da.on('failed', (reason) => {
      expect(reason).toBe('Failed to get authorisation info: Error: Invalid identifier "   "');
      done();
    });
    da.authorise('   ', ["9", "6", "10"]);
  });

  it('Authorisation failed - token is all whitespace', function(done) {
    var da = new DataAccess(constants.AUTH_TYPE_TOKEN);
    da.on('failed', (reason) => {
      expect(reason).toBe('Failed to get authorisation info: Error: Invalid identifier "   "');
      done();
    });
    da.authorise('   ', ["9", "6", "10"]);
  });

  it('Authorisation failed - email is expected', function(done) {
    config.provide(email_conf);
    var da = new DataAccess(constants.AUTH_TYPE_USER);
    da.on('failed', (reason) => {
      expect(reason).toBe('Failed to get authorisation info: Error: Invalid identifier "alice"');
      done();
    });
    da.authorise('alice', ["9", "6", "10"]);
  });

  it('Authorisation failed - email is incorrect', function(done) {
    config.provide(email_conf);
    var da = new DataAccess(constants.AUTH_TYPE_USER);
    da.on('failed', (reason) => {
      expect(reason).toBe('Failed to get authorisation info: Error: Invalid identifier "alice@boom.com"');
      done();
    });
    da.authorise('alice@boom.com', ["9", "6", "10"]);
  });

  it('Authorisation failed - access_control_group_id value is missing', function(done) {
    config.provide(email_conf);
    var da = new DataAccess(constants.AUTH_TYPE_USER);
    da.on('failed', (reason) => {
      expect(reason).toBe(
        'Failed to get authorisation info: Error: Some access group ids are not defined');
      done();
    });
    da.authorise('alice@boom.co.uk', ["9", "", "10"]);
  });

  it('Authorised - username', function(done) {
    config.provide(noemail_conf);
    var da = new DataAccess(constants.AUTH_TYPE_USER);
    da.on('authorised', () => {
       done();
    });
    da.on('failed', (reason) => {
      fail(reason);
      done();
    });
    da.authorise('alice', ['1', '2', '3']);
  });

  it('Authorised - token', function(done) {
    var da = new DataAccess(constants.AUTH_TYPE_TOKEN);
    da.on('authorised', () => {
      done();
    });
    da.on('failed', (reason) => {
      fail(reason);
      done();
    });
    da.authorise('abc', ['1', '2', '3']);
  });

  it('Authorisation failed - no auth for some of the files', function(done) {
    config.provide(noemail_conf);
    var da = new DataAccess(constants.AUTH_TYPE_USER);
    da.on('failed', (reason) => {
       expect(reason).toBe('Failed to get authorisation info: Error: Not authorised for those files');
       done();
    });
    da.authorise('alice', ["8", "6"]);
  });

  afterAll(function(done) {
    dummyAuthServ.close(done);
  });
});
