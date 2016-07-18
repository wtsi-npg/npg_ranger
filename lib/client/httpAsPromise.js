"use strict";

const assert    = require('assert');

const http      = require('http');
const url_tools = require('url');

const DEFAULT_APPLICATION_ERROR_CODE = 424;

/**
 * httpAsPromise module.
 * @module lib/client/httpAsPromise
 *
 * @author Jaime Tovar
 * @copyright Genome Research Limited 2016
 */
class HTTPAsPromise {

  /**
   * Create a HTTPAsPromise
   * @param {string} uri - uri of the resource
   * @param {object} [headers] - headers to pass in the request
   */
  constructor ( uri, headers ) {
    assert(uri, 'uri is required');
    if ( headers ) {
      assert.strictEqual(typeof headers, 'object', 'headers must be an object');
    }

    this.uri = uri;
    this.headers = headers || {};
    this.req = undefined;
  }

  _shorten ( theString, upTo ) {
    assert( theString, 'theString to shorten is required' );
    assert( typeof theString === 'string', 'theString must be of type string but came as: ' + typeof theString);
    upTo = upTo || 75;
    assert( typeof upTo === 'number', 'max length must be numeric' );

    if ( theString.length > upTo ) {
      theString = theString.substring(0, upTo - 3) + '...';
    } else {
      theString = theString.substring(0, upTo);
    }

    return theString;
  }

  _procEncoded( uri ) {
    assert(uri, 'uri is required');
    var regex = /^(data:[\S]+\/([\S]+);([\S]+),)(.*)/;

    var matches = uri.match(regex);
    var buffer;

    if ( matches ) {
      // var encString = matches[1];
      // var type      = matches[2];
      var encoding = matches[3];
      var data     = matches[4];

      buffer = new Buffer(data, encoding);
    } else {
      throw new Error('Unable to decode uri: ' + uri);
    }

    return buffer;
  }

  _reqOptionsForURL ( uri, agent ) {
    assert(uri, 'url is required');
    agent = agent || false;
    var parsedUrl = url_tools.parse( uri );

    var options = {
      agent: agent  // create a new agent just for this one request
    };
    if ( parsedUrl.hostname ) {
      options.hostname = parsedUrl.hostname;
    }
    if ( parsedUrl.port ) {
      options.port = parsedUrl.port;
    } else {
      options.port = parsedUrl.scheme === 'https' ? 443 : 80;
    }
    if ( parsedUrl.path ) {
      options.path = parsedUrl.path;
    }

    options.headers = this.headers;

    return options;
  }

  _buildReject( uri, status, message ) {
    /* console.log ('uri ' + uri +
                 ' status ' + status +
                 ' message ' + message); */
    let reason = {};
    reason.rejectUrl     = this._shorten( uri );
    reason.rejectStatus  = status;
    reason.rejectMessage = message;
    return reason;
  }

  /**
   * @return {Promise} Promise encapsulating the requests for data
   */
  run () {
    var self = this;

    var dataPromise = new Promise( ( resolve, reject ) => {
      if ( self.uri.startsWith('data:application') ) {
        self.req = {};
        let data;

        try {
          data = self._procEncoded( self.uri );
        } catch ( error ) {
          reject (
            self._buildReject(
              self.uri,
              DEFAULT_APPLICATION_ERROR_CODE,
              'Error while decoding data:application uri'
            )
          );
        }

        self.req.response      = data;
        self.req.headers       = { 'content-type': 'application/octet-stream' };
        self.req.rawHeaders    = [ 'content-type', 'application/octet-stream' ];
        self.req.trailers      = {};
        self.req.rawTrailers   = [];
        self.req.status        = 200;
        self.req.statusMessage = 'OK';
        self.req.url           = self.uri;
        self.req.readable      = true;

        resolve(self.req);
      } else {
        var options;
        try {
          options = self._reqOptionsForURL( self.uri );
        } catch (e) {
          reject (
            self._buildReject(
              self.uri,
              DEFAULT_APPLICATION_ERROR_CODE,
              e.toString()
            )
          );
        }

        var body = [];

        self.req = http.request(options, ( response ) => {
          response.on('data', ( data ) => {
            let dataBuffer = new Buffer(data, '');
            body.push(dataBuffer);
          });
          response.on('end', () => {
            // console.log('on end');
            let toCopy = 'headers rawHeaders trailers rawTrailers statusMessage'.split(' ');
            for ( let i = 0; i < toCopy.length; i++ ) {
              let name = toCopy[i];
              self.req[name] = response[name];
            }
            self.req.status = response.statusCode;
            self.req.response = Buffer.concat(body);
            self.req.url      = self.uri;
            self.req.readable = true;

            // TODO reject if trailer was set to truncated.

            if ( self.req.status === 200 || self.req.status === 206 ) {
              resolve(self.req);
            } else {
              reject (
                self._buildReject(
                  self.uri,
                  self.req.status,
                  self.req.statusMessage
                )
              );
            }
          });
        });

        // The browser implementation needs to know all returned data must be treated
        // as binary. Therefore I set responseType as 'arraybuffer' manually. The
        // default is binary for the node implementation.
        if ( self.req.xhr ) {
          self.req.xhr.responseType = 'arraybuffer';
        }

        self.req.on('error', ( error ) => {
          // console.log('on error');
          reject (
            self._buildReject(
              self.uri,
              DEFAULT_APPLICATION_ERROR_CODE,
              '' + error
            )
          );
        });

        // Create timeout for reject

        self.req.end();
      }
    });

    return dataPromise;
  }
}

module.exports = HTTPAsPromise;
