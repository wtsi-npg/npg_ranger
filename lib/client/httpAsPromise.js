"use strict";

const assert    = require('assert');

const http      = require('http');
const url_tools = require('url');

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

  _procEncoded( uri ) {
    var regex = /^(data:[\S]+\/([\S]+);([\S]+),)(.*)/;

    var matches   = uri.match(regex);
    // var encString = matches[1];
    // var type      = matches[2];
    var encoding  = matches[3];
    var data      = matches[4];
    var buffer    = new Buffer(data, encoding);

    return buffer;
  }

  _reqOptionsForURL ( uri, agent ) {
    assert(uri, 'url is required');
    agent = agent || false;
    var parsedUrl = url_tools.parse( uri );

    var options = {
      agent:    agent  // create a new agent just for this one request
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

  run () {
    var self = this;

    var dataPromise = new Promise( ( resolve, reject ) => {
      if ( self.uri.startsWith('data:application') ) {
        try {
          let data = self._procEncoded( self.uri );
          let req = {};

          req.headers       = { 'content-type': 'application/octet-stream' };
          req.rawHeaders    = [ 'content-type', 'application/octet-stream' ];
          req.trailers      = {};
          req.rawTrailers   = [];
          req.statusCode    = 200;
          req.statusMessage = 'OK';
          req.response      = data;
          req.url           = self.uri;
          req.readable      = true;

          resolve(req);
        } catch ( error ) {
          reject(error);
        }
      } else {
        var options = self._reqOptionsForURL( self.uri );
        var body = [];

        self.req = http.request(options, ( response ) => {
          response.on('data', ( data ) => {
            let dataBuffer = new Buffer(data, '');
            body.push(dataBuffer);
          });
          response.on('end', () => {
            console.log('End ' + self.uri);
            let toCopy = 'headers rawHeaders trailers rawTrailers statusCode statusMessage'.split(' ');
            for ( let i = 0; i < toCopy.length; i++ ) {
              let name = toCopy[i];
              self.req[name] = response[name];
            }
            self.req.response      = Buffer.concat(body);
            self.req.url           = self.uri;
            self.req.readable      = true;
            console.log('Got ' + self.req.response.byteLength + ' bytes');
            if ( response.statusCode === 200 || response.statusCode === 206 ) {
              resolve(self.req);
            } else {
              reject(self.req);
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
          reject(error);
        });

        self.req.end();
      }
    });

    return dataPromise;
  }
}

module.exports = HTTPAsPromise;
