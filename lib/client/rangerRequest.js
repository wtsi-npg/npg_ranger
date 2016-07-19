"use strict";

const EventEmitter  = require('events');
const assert        = require('assert');

const LOGGER        = require('winston');

const HTTPAsPromise = require('./httpAsPromise');

const UNSENT          = 0; // Client has been created. open() not called yet.
const OPENED          = 1; // open() has been called.
// const HEADERS_RECEIVED = 2; // send() has been called, and headers and status are available.
// const LOADING          = 3; // Downloading; responseText holds partial data.
const DONE            = 4; // The operation is complete.

const ENC_DATA_PREFIX = 'data:application/vnd-ga4gh;base64,';

const DEFAULT_APPLICATION_ERROR_CODE = 424;

/**
 * npgRanger module.
 * @module lib/client/npgRanger
 *
 * @author Jaime Tovar
 * @copyright Genome Research Limited 2016
 */
class RangerRequest extends EventEmitter {

  /**
   * Create a RangerRequest
   */
  constructor () {
    super();

    this.readyState       = UNSENT;
    this.status           = 0;
    this.response         = '';
    this._responseHeaders = {};
    this.withCredentials  = false;
    this.requestHeaders   = {};
    this._requestOptions  = {};

    // TODO add max number or retries
    // TODO add max number of redirects

    var self = this;
    this.on('readystatechange', () => { self.onreadystatechange(); });
    this.on('error',            () => { self.onerror(); });
  }

  /**
   * Set up the configuration for the request
   * @param {string} method - The http method to use for the request
   * @param {string} url - The url of the resource
   * @param {boolean} [async] - The request will be async, by default true
   */
  open( method, url, async ) {
    assert(method, 'method is required');
    assert(url, 'url is required');
    assert(typeof method === 'string', 'method must be a string');
    assert(typeof url === 'string', 'url must be a string');

    if ( typeof async !== 'undefined' ) {
      assert(typeof async === 'boolean', 'async can only be boolean');
      assert(async === true, 'Only async requests are supported');
    } else {
      async = true;
    }

    this._requestOptions.url    = url;
    this._requestOptions.method = method;
    this._requestOptions.async  = async;

    this._setReadyState(OPENED);
  }

  /*
   * Merge properties of two objects passed as parameters, give priority
   * (don't overwrite) properties of first parameter.
   * Returns a new object with the merged properties
   * @param {object} h1 - First header object, properties of this object will
   *                      have priority over the second object
   * @param {object} h2 - Second header object
   * @return {object} New object with merged properties
   */
  _mergeHeaders(h1, h2) {
    assert(h1, 'header 1 is required');
    assert(h2, 'header 2 is required');
    var mergedHeaders = {};

    for ( let i in h1 ) {
      mergedHeaders[i] = h1[i];
    }
    for ( let i in h2 ) {
      if ( !mergedHeaders[i] ) {
        mergedHeaders[i] = h2[i];
      }
    }
    return mergedHeaders;
  }

  send( payload ) {
    payload = payload || '';
    var self = this;

    if ( this.readyState !== OPENED ) {
      throw new Error('The object state must be OPENED');
    }

    // if ( this.withCredentials ) {
    //   this.requestHeaders.withCredentials = this.withCredentials;
    // }

    LOGGER.debug(self._requestOptions.url);
    var httpP = new HTTPAsPromise( self._requestOptions.url, this.requestHeaders );
    var dataPromise = httpP.run();

    dataPromise.then( ( r ) => {
      var contentType = self._getContentType(r);
      if ( contentType.startsWith('application/json') ) {
        try {
          self._procJSON(r.response);
        } catch (e) {
          LOGGER.error(e);
          self.status        = DEFAULT_APPLICATION_ERROR_CODE;
          self.statusMessage = "Error while processing JSON " + e;
          delete self.response;
          self._setReadyState(DONE);
        }
      } else {
        LOGGER.debug('Data in response');
        self.status        = r.status;
        self.statusMessage = r.statusMessage;
        self.response      = r.response;
        self._setReadyState(DONE);
      }
    }, ( reason ) => {
      self.status        = reason.rejectStatus;
      self.statusMessage = reason.rejectMessage + ' for ' + reason.rejectUrl;
      delete self.response;
      self._setReadyState(DONE);
    });
  }

  _getContentType( r ) {
    var contentType = r.headers['content-type'];
    contentType = contentType && contentType.toLowerCase ? contentType.toLowerCase() : '';
    return contentType;
  }

  _procJSON( data ) {
    var decdStr;
    decdStr = String.fromCharCode.apply(null, new Uint8Array( data ));

    var jsonResponse = JSON.parse(decdStr);
    LOGGER.debug('Response with JSON: ' + decdStr);
    var uris         = [];
    var headers4uris = [];
    var uriPromises  = [];

    var prefix;
    var suffix;

    for (let i = 0; i < jsonResponse.urls.length; i++ ) {
      let uri;
      let headers = this.requestHeaders;
      if ( jsonResponse.urls[i].url ) { // TODO use only this one once all servers are updated.
        uri = jsonResponse.urls[i].url;
        LOGGER.debug('Found URI from JSON ' + uri);
        if ( jsonResponse.urls[i].headers ) {

          // Merge headers, give priority to new headers provided
          headers = this._mergeHeaders(jsonResponse.urls[i].headers, headers);
        }
      } else {
        // uri = jsonResponse.urls[i];
        // TODO Replace as soon as CH updates his server
        throw new Error('Malformed JSON redirect, missing url field in: ' + JSON.stringify(jsonResponse.urls[i]));
      }
      uris[i] = uri;
      headers4uris[i] = headers;
    }

    if ( jsonResponse.prefix && jsonResponse.prefix.length ) {
      prefix = ENC_DATA_PREFIX + jsonResponse.prefix;
      uris.unshift(prefix);
    }

    if ( jsonResponse.suffix && jsonResponse.suffix.length ) {
      suffix = ENC_DATA_PREFIX + jsonResponse.suffix;
      uris.push(suffix);
    }

    for ( let i = 0; i < uris.length; i++ ) {
      let httpP = new HTTPAsPromise( uris[i], headers4uris[i] );
      uriPromises[i] = httpP.run();
    }

    Promise.all(uriPromises).then( values => {
      try {
        let dataChunks = [];
        for ( let i = 0; i < values.length; i ++ ) {
          dataChunks[i] = values[i].response;
        }

        let buffer = Buffer.concat(dataChunks);
        // Browser implementation of zlib is expecting an ArrayBuffer. The node version
        // can directly work with a Buffer object.
        if ( buffer.toArrayBuffer ) {
          buffer = buffer.toArrayBuffer();
        }
        this.response      = buffer;
        this.status        = 200;
        this.statusMessage = 'OK';
      } catch (e) {
        this.response      = 'Error while building reponse ' + e;
        this.status        = DEFAULT_APPLICATION_ERROR_CODE;
        this.statusMessage = 'Error while building reponse ' + e;
      }
      this._setReadyState(DONE);
    }, ( reason ) => {
      let url           = reason.rejectUrl;
      let status        = reason.rejectStatus;
      let statusMessage = reason.rejectMessage;
      this.response = 'Error while requesting resources provided in JSON ' +
                      url + ' ' +
                      status + ' ' +
                      statusMessage;
      this.response.status = DEFAULT_APPLICATION_ERROR_CODE;
      this._setReadyState(DONE);
    });
  }

  _cloneResponse(response) {
    var newResponse = {};
    var toCopy = 'headers rawHeaders rawTrailers trailers readable statusCode statusMessage url'.split(' ');
    for ( let i = 0; i < toCopy.length; i++ ) {
      newResponse[toCopy[i]] = response[toCopy[i]];
    }
    return newResponse;
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

  getResponseHeader( headerName ) {
    assert(headerName, 'headerName is required');
    assert(headerName.toLowerCase, 'headerName must be a string');
    headerName = headerName.toLowerCase();
    return this._responseHeaders[headerName];
  }

  onreadystatechange() {}
  onerror() {}
}

module.exports = RangerRequest;
