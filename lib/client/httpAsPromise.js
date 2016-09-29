"use strict";

const assert    = require('assert');
const http      = require('http');
const https     = require('https');
const url_tools = require('url');
// const LOGGER      = require('winston');

const DEFAULT_APPLICATION_ERROR_CODE = 424;

/**
 * The build in Promise
 * @external Promise
 * @see  {@link https://developer.mozilla.org/en/docs/Web/JavaScript/Reference/Global_Objects/Promise|Promise}
 */

/**
 * NodeJS implementation of assert
 * @external assert
 * @see  {@link https://nodejs.org/dist/latest-v4.x/docs/api/assert.html|assert}
 */

/**
 * NodeJS implementation of http
 * @external http
 * @see  {@link https://nodejs.org/dist/latest-v4.x/docs/api/http.html|http}
 */

/**
 * NodeJS implementation of https
 * @external https
 * @see  {@link https://nodejs.org/dist/latest-v4.x/docs/api/https.html|https}
 */

/**
 * NodeJS implementation of url
 * @external url
 * @see  {@link https://nodejs.org/dist/latest-v4.x/docs/api/url.html|url}
 */

/**
 * httpAsPromise module.
 * @module client/httpAsPromise
 *
 * @requires {@link external:assert|assert}
 * @requires {@link external:http|http}
 * @requires {@link external:https|https}
 * @requires {@link external:url|url}
 *
 * @author Jaime Tovar
 * @copyright Genome Research Limited 2016
 */

/**
 * A {@link external:Promise|Promise} encapsulating an http request.
 *
 * @example
 * var HTTPAsPromise  = require('./httpAsPromise');
 * var url            = 'http://192.168.0.1:5050/' +
 *                      'directresources/AA0011?r=1&range=165000-175000&format=BAM';
 *
 * var httpP = new HTTPAsPromise( url );
 * var dataPromise = httpP.run();
 * dataPromise.then( ( r ) => {
 *   var contentType = r.headers['content-type'];
 *   contentType = contentType && contentType.toLowerCase ? contentType.toLowerCase() : '';
 *   if ( contentType.startsWith('application/json') ) {
 *     try {
 *       console.log(r.response);
 *     } catch (e) {
 *       console.log(e);
 *     }
 *   } else {
 *     console.log('Unexpected content type');
 *   }
 * }, ( reason ) => {
 *   console.log( reason.rejectStatus + ' ' +
 *                reason.rejectMessage + ' for ' +
 *                reason.rejectUrl );
 * });
 *
 */
class HTTPAsPromise {

  /**
   * @param {string} uri - uri of the resource
   * @param {object} [headers] - headers to pass as part of the request. Defaults to
   *                             and empty object
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
      options.port = parsedUrl.protocol === 'https:' ? 443 : 80;
    }
    if ( parsedUrl.protocol === 'https:' ) {
      options.protocol = parsedUrl.protocol;
    }
    if ( parsedUrl.path ) {
      options.path = parsedUrl.path;
    }

    options.headers = this.headers;
    // http-browserify will set it as true by default just before the request
    // therefore I make sure undefined and false go to false.
    if ( !options.withCredentials ) {
      options.withCredentials = false;
    }

    return options;
  }

  _buildReject( uri, status, message ) {
    // LOGGER.debug('uri ' + uri +
    //             ' status ' + status +
    //             ' message ' + message);
    let reason = {};
    reason.rejectUrl     = this._shorten( uri );
    reason.rejectStatus  = status;
    reason.rejectMessage = message;
    return reason;
  }

  /**
   * Run the request inside a promise
   * @return {@link external:Promise|Promise} encapsulating the
   *         requests for data
   */
  run () {
    var self = this;

    var dataPromise = new Promise( ( resolve, reject ) => {
      if ( self.uri.startsWith('data:') ) {
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

        let requestProtocol = options.protocol === 'https:' ? https : http;
        self.req = requestProtocol.request(options, ( response ) => {
          response.on('data', ( data ) => {
            // LOGGER.debug('on data');
            let dataBuffer = new Buffer(data, '');
            body.push(dataBuffer);
          });
          response.on('end', () => {
            // LOGGER.debug('on end');
            let toCopy = 'headers rawHeaders trailers rawTrailers statusMessage'.split(' ');
            for ( let i = 0; i < toCopy.length; i++ ) {
              let name = toCopy[i];
              self.req[name] = response[name];
            }
            self.req.status   = response.statusCode;
            self.req.response = Buffer.concat(body);
            self.req.url      = self.uri;
            self.req.readable = true;

            if ( response.trailers['data-truncated'] && response.trailers['data-truncated'] === 'true' ) {
              reject(
                self._buildReject(
                  self.uri,
                  self.req.status,
                  'Data was truncated'
                  )
                );
            } else if ( self.req.status === 200 || self.req.status === 206 ) {
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
          // LOGGER.error('on error ' + error)
          console.log('on error');
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

/**
 * HTTPAsPromise An http request encapsulated as {@link external:Promise|Promise}
 * @type {HTTPAsPromise}
 */
module.exports = HTTPAsPromise;
