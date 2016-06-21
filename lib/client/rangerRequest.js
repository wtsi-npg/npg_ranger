"use strict";

const EventEmitter = require('events');
const assert       = require('assert');

// For http
const http      = require('http');
const url_tools = require('url');

const UNSENT           = 0; // Client has been created. open() not called yet.
const OPENED           = 1; // open() has been called.
// const HEADERS_RECEIVED = 2; // send() has been called, and headers and status are available.
// const LOADING          = 3; // Downloading; responseText holds partial data.
const DONE             = 4; // The operation is complete.

class RangerRequest extends EventEmitter {

  constructor () {
    super();

    this.readyState      = UNSENT;
    this.status          = 0;
    this.reponseType     = '';
    this.response        = '';
    this.withCredentials = false;
    this.requestHeaders  = {};
    this._requestOptions = {};

    var self = this;
    this.on('readystatechange', () => { self.onreadystatechange(); });
    this.on('error',            () => { self.onerror(); });
  }

  open( method, url, async ) {
    assert(url, 'url is required');
    assert(method, 'method is required');
    async = async || true;

    this._requestOptions.url    = url;
    this._requestOptions.method = method;
    this._requestOptions.async  = async;

    this._setReadyState(OPENED);
  }

  _reqOptionsForURL ( url, agent ) {
    assert(url, 'url is required');
    agent = agent || false;
    var parsedUrl = url_tools.parse( url );

    var options = {
      path:     parsedUrl.path,
      hostname: parsedUrl.hostname,
      port:     parsedUrl.port,
      agent:    agent,  // create a new agent just for this one request
    };
    return options;
  }

  _reqOptionsForSocket ( socketPath, path, user ) {
    assert( socketPath, 'socketPath is required' );
    assert( path, 'path is required');
    var options = {
      socketPath: socketPath,
      path:       path
    };
    if ( user ) {
      options.headers = { 'X-Remote-User': user };
    }
    return options;
  }

  send( payload ) {
    payload = payload || '';
    var self = this;

    if ( this.readyState !== OPENED ) {
      throw new Error('The object state must be OPENED');
    }

    var req;

    var dataPromise = new Promise( ( resolve, reject ) => {
      var body    = '';
      var options = self._reqOptionsForURL( self._requestOptions.url );

      req = http.get(options , ( response ) => {

        // Continuously update stream with data
        response.on('data', ( data ) => {
          body += data;
        });
        response.on('end', () => {

          // Data reception is done, do whatever with it!
          // console.log('RESPONSE CODE ' + response.statusCode + ' ' + response.statusMessage);
          // console.log('RAW HEADERS ' + response.rawHeaders);
          // console.log('RAW TRAILERS ' + response.rawTrailers);
          // console.log('TRAILER data-truncated' +
          //  (JSON.stringify(response.trailers)));
          // console.log("RECEIVED: " + body);

          if ( response.statusCode === 200 || response.statusCode === 206 ) {
            self.response = body;
            resolve();
          } else {
            self.response = response.statusMessage;
            self.response.statusCode = response.statusCode;
            reject();
          }
        });
      });
    });

    req.on('error', (error) => {
      self.response = error.message;
      self.response.statusCode = error.status;
      self._setReadyState(DONE);
      dataPromise.reject(error);
    });

    dataPromise.then( () => {
      self._setReadyState(DONE);
    });
  }

  _setReadyState( readyState ) {
    assert(readyState, 'readyState is required');
    this.readyState = readyState;
    this.emit('readystatechange');
  }

  setRequestHeader( name, value ) {
    assert(name, 'name is required');
    assert(value, 'value is required');
    this.requestHeaders[name] = value;
  }

  getResponseHeader() {
    return '';
  }

  onreadystatechange() {}
  onerror() {}
}

module.exports = RangerRequest;
