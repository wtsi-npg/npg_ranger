/* globals describe, it, expect, beforeAll, afterAll, beforeEach, afterEach */

'use strict';

const https       = require('https');
const url         = require('url');
const fs          = require('fs');

const DataAccess  = require('../../lib/server/auth.js');
const config      = require('../../lib/config.js');
const constants   = require('../../lib/constants.js');
const test_utils  = require('./test_utils');

let BASE_PORT  = 14000;
let PORT_RANGE = 200;
let PORT       = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;

let tmpDir = config.tempFilePath('npg_ranger_auth_test_');

let _handle_request = ( req, res ) => {
  let auth_groups = ['1', '2', '3'];
  let test_groups = groups => {
    return groups.every( individual_group => {
      return individual_group.some( g_element => {
         return auth_groups.indexOf(g_element) != -1;
      });
    });
  };

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
    console.log(JSON.stringify(reqdata));
    if (requrl.pathname === constants.AUTH_URL_TOKEN) {
      ok = (reqdata.token === 'abc' && test_groups(reqdata.groups));
      res.end(JSON.stringify({ok: ok}));
    } else if (requrl.pathname === constants.AUTH_URL_USER) {
      ok = (reqdata.user === 'alice' && test_groups(reqdata.groups));
      res.end(JSON.stringify({ok: ok}));
    } else {
      res.statusCode = 404;
      res.end();
    }
  });
};

let create_https_server = (cert, key, ca) => {
  let options = {
    key:  fs.readFileSync(key),
    cert: fs.readFileSync(cert),
    ca:   fs.readFileSync(ca)
  };
  let server = https.createServer(options, _handle_request);
  server.listen(PORT);
  return server;
};

['https'].forEach( serverType => {
  describe(`Authorisation with auth server on ${serverType}`, function() {
    let closeServ;
    let dummyAuthServ;

    let serv_cert;
    let serv_key;
    let ca;
    let client_cert;
    let client_key;

    let addAuthSSL = opts => {
      opts.auth_cert = client_cert;
      opts.auth_key  = client_key;
      opts.auth_ca   = ca;

      return opts;
    };

    let noemail_conf = function() {
      let opts = {
        emaildomain: null,
        authurl: `${serverType}://localhost:${PORT}/`
      };
      return addAuthSSL(opts);
    };

    let email_conf = function() {
      let opts = {
        emaildomain: 'boom.co.uk',
        authurl: `${serverType}://localhost:${PORT}/`
      };
      return addAuthSSL(opts);
    };

    beforeAll( (done) => {
      // Provide a minimal standin for npg_sentry.
      // Authorises user 'alice' or token 'abc' for groups 1, 2 and 3.
      let ca_prefix    = 'ca';
      let cert1_prefix = 'serv';
      let cert2_prefix = 'client';
      test_utils.create_certificates(tmpDir, ca_prefix, cert1_prefix, cert2_prefix, () => {
        ca          = `${tmpDir}/${ca_prefix}.cert`;
        serv_key    = `${tmpDir}/${cert1_prefix}.key`;
        serv_cert   = `${tmpDir}/${cert1_prefix}.cert`;
        client_key  = `${tmpDir}/${cert2_prefix}.key`;
        client_cert = `${tmpDir}/${cert2_prefix}.cert`;
        config.provide(noemail_conf);
        dummyAuthServ = create_https_server(serv_cert, serv_key, ca);
        closeServ = dummyAuthServ.close;
        done();
      });
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
        expect(reason).toMatch(/Unexpected error while processing authorisation/i);
        done();
      });
      da.authorise('   ', [["9"], ["6"], ["10"]]);
    });

    it('Authorisation failed - token is all whitespace', function(done) {
      var da = new DataAccess(constants.AUTH_TYPE_TOKEN);
      da.on('failed', (reason) => {
        expect(reason).toMatch(/Unexpected error while processing authorisation/i);
        done();
      });
      da.authorise('   ', [["9"], ["6"], ["10"]]);
    });

    it('Authorisation failed - email is expected', function(done) {
      config.provide(email_conf);
      var da = new DataAccess(constants.AUTH_TYPE_USER);
      da.on('failed', (reason) => {
        expect(reason).toMatch(/Unexpected error while processing authorisation/i);
        done();
      });
      da.authorise('alice', [["9"], ["6"], ["10"]]);
    });

    it('Authorisation failed - email is incorrect', function(done) {
      config.provide(email_conf);
      var da = new DataAccess(constants.AUTH_TYPE_USER);
      da.on('failed', (reason) => {
        expect(reason).toMatch(/Unexpected error while processing authorisation/i);
        done();
      });
      da.authorise('alice@boom.com', [["9"], ["6"], ["10"]]);
    });

    [
      [["9"], [""], ["10"]],
      [["9"], [], ["10"]]
    ].forEach( groups => {
      it('Authorisation failed - access_control_group_id value is missing', function(done) {
        config.provide(email_conf);
        var da = new DataAccess(constants.AUTH_TYPE_USER);
        da.on('failed', (reason) => {
          expect(reason).toMatch(/Not authorised for those files/i);
          done();
        });
        da.authorise('alice@boom.co.uk', groups);
      });
    });

    [
      ['1', '2', '3'],
      [['1'], ['2'], ['3']],
      [['1'], '2', ['3']]
    ].forEach( access_groups => {
      it('Authorised - username array groups', function(done) {
        config.provide(noemail_conf);
        var da = new DataAccess(constants.AUTH_TYPE_USER);
        da.on('authorised', () => {
           done();
        });
        da.on('failed', (reason) => {
          console.log(reason);
          done.fail(reason);
        });
        da.authorise('alice', access_groups);
      });
    });

    [
      ['1', '2', '3'],
      [['1'], ['2'], ['3']],
      [['1'], '2', ['3']]
    ].forEach( access_groups => {
      it('Authorised - token array groups', function(done) {
        var da = new DataAccess(constants.AUTH_TYPE_TOKEN);
        da.on('authorised', () => {
          done();
        });
        da.on('failed', (reason) => {
          done.fail(reason);
        });
        da.authorise('abc', access_groups);
      });
    });

    it('Authorisation failed - no auth for any of the files', function(done) {
      config.provide(noemail_conf);
      var da = new DataAccess(constants.AUTH_TYPE_USER);
      da.on('failed', (reason) => {
        expect(reason).toMatch(/Not authorised for those files/i);
        done();
      });
      da.authorise('alice', [["8"], ["6"]]);
    });

    [
      [["1"], ["6"]],
      ["1", ["6"]],
      ["1", ["2"], ["6"]],
      ["1", ["2"], "6"],
    ].forEach( access_groups => {
      it('Authorisation failed - no auth for some of the files', function(done) {
        config.provide(noemail_conf);
        var da = new DataAccess(constants.AUTH_TYPE_USER);
        da.on('failed', (reason) => {
          expect(reason).toMatch(/Not authorised for those files/i);
          done();
        });
        da.authorise('alice', access_groups);
      });
    });

    afterAll(function(done) {
      dummyAuthServ.close(done);
    });
  });
});

describe('Validate extra properties in response', () => {
  let server;

  let serv_cert;
  let serv_key;
  let ca;
  let client_cert;
  let client_key;

  let BAD_SERVER_PORT = Math.floor(Math.random() * PORT_RANGE) + BASE_PORT;

  let noemail_conf = function() {
    let opts = {
      emaildomain: null,
      authurl: `https://localhost:${BAD_SERVER_PORT}/`
    };
    return addAuthSSL(opts);
  };

  let addAuthSSL = opts => {
    opts.auth_cert = client_cert;
    opts.auth_key  = client_key;
    opts.auth_ca   = ca;

    return opts;
  };

  let _bad_handle_request = ( req, res ) => {
    req.on('data', () => {});
    req.on('end', () => {
      res.statusCode = 200;
      res.end(JSON.stringify({ok: true, notOK: false}));
    });
  };

  let createBadPropsServer = (cert, key, ca) => {
    let options = {
      key:  fs.readFileSync(key),
      cert: fs.readFileSync(cert),
      ca:   fs.readFileSync(ca)
    };
    server = https.createServer(options, _bad_handle_request);
    server.listen(BAD_SERVER_PORT);
  };

  beforeEach( (done) => {
    let ca_prefix    = 'ca';
    let cert1_prefix = 'serv';
    let cert2_prefix = 'client';
    test_utils.create_certificates(tmpDir, ca_prefix, cert1_prefix, cert2_prefix, () => {
      ca          = `${tmpDir}/${ca_prefix}.cert`;
      serv_key    = `${tmpDir}/${cert1_prefix}.key`;
      serv_cert   = `${tmpDir}/${cert1_prefix}.cert`;
      client_key  = `${tmpDir}/${cert2_prefix}.key`;
      client_cert = `${tmpDir}/${cert2_prefix}.cert`;
      config.provide(noemail_conf);
      createBadPropsServer(serv_cert, serv_key, ca);
      done();
    });
  });

  afterEach( (done) => {
    server.close();
    done();
  });

  [['1', '2', '3'], [['1'], ['2'], ['3']]].forEach( group_set  => {
    it('Fails with extra properties in response', (done) => {
      config.provide(noemail_conf);
      var da = new DataAccess(constants.AUTH_TYPE_USER);
      da.on('authorised', () => {
        console.log("HERE");
        done.fail('unexpected authorised');
      });
      da.on('failed', (reason) => {
        expect(reason).toMatch(/Unexpected error while processing authorisation/i);
        done();
      });
      da.authorise('alice', group_set);
    });
  });
});
