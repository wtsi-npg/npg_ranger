"use strict";

const EventEmitter = require('events');
const assert       = require('assert');

// const LOGGER    = require('winston');

const HTTPAsPromise = require('./httpAsPromise');

/**
 * @external events
 * @see {@link https://nodejs.org/dist/latest-v4.x/docs/api/events.html|events}
 */

/**
 * @external EventEmitter
 * @see {@link https://nodejs.org/dist/latest-v4.x/docs/api/events.html#events_class_eventemitter|EventEmitter}
 */

/**
 * <p>Module which provides the HTTPAsPromise object and related utils.</p>
 *
 * <p>Defines implemented ready states reported by HTTPAsPromise object during a
 * request
 * (see <a href="https://developer.mozilla.org/en-US/docs/Web/API/XMLHttpRequest/readyState">
 * XMLHttpRequest.readyState</a>).</p>
 *
 * <p>Defines http error code used to report npg_ranger library-exclusive
 * errors (<code>424</code>)</p>
 *
 * @module client/rangerRequest
 *
 * @requires {@link external:events|events}
 * @requires {@link external:assert|assert}
 * @requires module:client/httpAsPromise
 *
 * @author Jaime Tovar
 * @copyright Genome Research Limited 2016
 */

/**
 * Set to <i>0</i>. Client has been created. <code>open()</code> not called yet.
 * @type {Number}
 */
const UNSENT = 0;

/**
 * Set to <i>1</i>. Matching XMLHttpRequest.<code>open()</code> has been called for this request.
 * @type {Number}
 */
const OPENED = 1;

// const HEADERS_RECEIVED = 2; // send() has been called, and headers and status are available.
// const LOADING          = 3; // Downloading; responseText holds partial data.

/**
 * Set to <i>4</i>. The operation is complete
 * @type {Number}
 */
const DONE   = 4;

const ENC_DATA_PREFIX = 'data:application/vnd-ga4gh;base64,';

/**
 * The 3-digit HTTP request error code. Set as <code>424</code> following the
 * 400's for client requests.
 * <p>This is the kind of error rised when for example:</p>
 * <ul>
 *  <li>
 *    there was an error but no error code can be obtained from the lower
 *    level http request implementation
 *  </li>
 *  <li>
 *    the process fails before the lower level http request is created
 *  </li>
 *  <li>
 *    the process fails after the lower level http request is finished but
 *    before this module can process it
 *  </li>
 * </ul>
 * @type {Number}
 */
const DEFAULT_APPLICATION_ERROR_CODE = 424;

/**
 * Class wrapping the business logic of a GA4GH request and exposing it, as
 * much as possible, as a {@link external:EventEmitter|EventEmitter} based on
 * the http request pattern.
 *
 * @example
 * var RangerRequest = require('npg_ranger').RangerRequest;
 * var url = 'http://192.168.0.1:5050/' +
 *           'resources/AA0011?referenceName=1&start=165000&end=175000&format=BAM';
 * var request = new RangerRequest();
 * request.open('GET', url);
 * // Add listeners
 * request.onreadystatechange = () => {
 *   console.log('readystatechange' + request.readyState);
 *   // ...
 * }
 * request.send('');
 *
 * @extends {external:EventEmitter}
 */
class RangerRequest extends EventEmitter {
  /**
   * Returns a RangerRequest object initilised with readyState set to
   * <code>UNSENT</code>. It also sets up default calls for error and
   * readystatechange events. It is expected listeners should will assigned
   * to process these events.
   */
  constructor () {
    super();

    this.readyState       = UNSENT;
    this.status           = 0;
    this.response         = '';
    this.withCredentials  = false;
    this.requestHeaders   = {};
    this._requestOptions  = {};
    this._responseHeaders = {};

    // TODO add max number or retries
    // TODO add max number of redirects

    var self = this;

    // Set default empty listeners for basic functionality
    this.on('readystatechange', () => { self.onreadystatechange(); });
    this.on('error',            () => { self.onerror(); });

    // Not yet emitted but will.
    this.on('abort',            () => { self.onabort(); });
    this.on('load',             () => { self.onload(); });
    this.on('loadstart',        () => { self.onloadstart(); });
    this.on('progress',         () => { self.onprogress(); });
  }

  /**
   * <p>Set up the configuration for the request. <code>method</code> and <code>url</code>
   * are passed directly to the internal implementation of the request. This object only
   * enforces the async parameter to be <code>true</code> (and is the default value). Sync
   * requests are not appropiate for the use case.</p>
   *
   * <p>Method is not enforced at this level. But only <code>'GET'</code> has been tested.</p>
   *
   * @example
   * // ...
   * var request = new RangerRequest();
   * request.open('GET', url, true);
   * request.send('');
   * // ...
   * var request2 = new RangerRequest();
   * request2.open('GET', url);
   * request2.send('');
   * // ...
   *
   * @param {string} method - The http method to use for the request. Usually a
   *                          <code>GET</code> is expected but not enforced.
   * @param {string} url - The url of the resource
   * @param {boolean} [async] - By default the request will be async. If a value is
   *                            passed, the value must be <code>true</code>.
   *                            Otherwise an <code>AssertionError</code> will be rised.
   *                            This parameter is part of the method's
   *                            signature just to maintain a signature compatible with
   *                            the general <code>request</code> object.
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
    assert(h1, 'headers 1 is required');
    assert(h2, 'headers 2 is required');
    assert(typeof h1 === 'object', 'headers 1 must be an object');
    assert(typeof h2 === 'object', 'headers 2 must be an object');
    var mergedHeaders = {};

    for ( let i in h1 ) {
      mergedHeaders[i] = h1[i];
    }
    for ( let i in h2 ) {
      if ( typeof mergedHeaders[i] === 'undefined' ) {
        mergedHeaders[i] = h2[i];
      }
    }
    return mergedHeaders;
  }

  /**
   * Sends the request with the specified payload. Usually an empty payload is expected.
   *
   * @example
   * // ...
   * request.open('GET', url);
   * request.send('');
   * // ...
   * request2.open('GET', url);
   * request2.send(null);
   * // ...
   * request3.open('GET', url);
   * request3.send();
   *
   * @param {object} [payload] - data to be sent as the body of the request. Likely
   *                             to be empty () or empty string (''). defaults to
   *                             empty string.
   */
  send( payload ) {
    payload = payload || '';
    var self = this;

    if ( this.readyState !== OPENED ) {
      throw new Error('The object state must be OPENED');
    }

    if ( this.withCredentials ) {
      this.requestHeaders.withCredentials = this.withCredentials;
    }

    // LOGGER.debug(self._requestOptions.url);
    var httpP = new HTTPAsPromise( self._requestOptions.url, this.requestHeaders );
    var dataPromise = httpP.run();

    dataPromise.then( ( r ) => {
      var contentType = self._getContentType(r);
      if ( contentType.startsWith('application/json') ) {
        try {
          self._procJSON(r.response);
        } catch (e) {
          // LOGGER.error(e);
          console.log(e);
          self.status        = DEFAULT_APPLICATION_ERROR_CODE;
          self.statusMessage = "Error while processing JSON " + e;
          delete self.response;
          self._setReadyState(DONE);
        }
      } else {
        // LOGGER.debug('Data in response');
        self.status        = r.status;
        self.statusMessage = r.statusMessage;
        let response = r.response;
        // Browser implementation of zlib is expecting an ArrayBuffer. The node version
        // can directly work with a Buffer object.
        if ( response.toArrayBuffer ) {
          response = response.toArrayBuffer();
        }
        self.response = response;
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
    // LOGGER.debug('Response with JSON: ' + decdStr);
    var uris         = [];
    var headers4uris = [];
    var uriPromises  = [];

    var prefix;
    var suffix;

    for (let i = 0; i < jsonResponse.urls.length; i++ ) {
      let uri;
      let headers = this.requestHeaders;
      if ( jsonResponse.urls[i].url ) {
        uri = jsonResponse.urls[i].url;
        // LOGGER.debug('Found URI from JSON ' + uri);
        if ( jsonResponse.urls[i].headers ) {

          // Merge headers, give priority to new headers provided
          headers = this._mergeHeaders(jsonResponse.urls[i].headers, headers);
        }
      } else {
        throw new Error( 'Malformed JSON redirect, missing url field in: ' +
                         JSON.stringify(jsonResponse.urls[i]) );
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
      this.status       = DEFAULT_APPLICATION_ERROR_CODE;
      let statusMessage = 'Error while requesting resources provided in JSON ' +
                          url + ' ' +
                          reason.rejectStatus + ' ' +
                          reason.rejectMessage;
      this.response      = statusMessage;
      this.statusMessage = statusMessage;
      this._setReadyState(DONE);
    });
  }

  _cloneResponse( response ) {
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

  /**
   * Set or overwrite the header for the request
   * @param {string} name - name of the header
   * @param {object} value - new value for the header
   */
  setRequestHeader( name, value ) {
    assert(name, 'name is required');
    assert(value, 'value is required');
    this.requestHeaders[name] = value;
  }

  /**
   * Get a header from the reponse by header name
   * @param  {string} headerName - name of the header required
   * @return {string} value of the header if it exists undefined otherwise
   */
  getResponseHeader( headerName ) {
    assert(headerName, 'headerName is required');
    assert(headerName.toLowerCase, 'headerName must be a string');
    headerName = headerName.toLowerCase();
    return this._responseHeaders[headerName];
  }

  /**
   * Empty method, expected to be overwritten by user. Will be called when
   * object raises a <code>readystatechange</code> event.
   */
  onreadystatechange() {}

  /**
   * Empty method, expected to be overwritten by user. Will be called when
   * object raises a <code>error</code> event.
   */
  onerror() { console.log('Default on.error call'); }

  /**
   * Waiting for implementation
   */
  onabort() { console.log('Default on.abort call'); }
  /**
   * Waiting for implementation
   */
  onload() { console.log('Default on.load call'); }
  /**
   * Waiting for implementation
   */
  onloadstart() { console.log('Default on.loadstart call'); }
  /**
   * Waiting for implementation
   */
  onprogress() { console.log('Default on.progress call'); }
}

module.exports = RangerRequest;
