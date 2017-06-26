"use strict";

/**
 * controller module.
 * @module server/controller
 *
 * @description Controller module for the streaming HTTP server.
 *
 * @example <caption>Example usage of the controller.</caption>
 *   const RangerController = require('../lib/server/controller.js');
 *   // Create a new object, no authorisation.
 *   let c = new RangerController(request, response, db);
 *   // handle request
 *   c.handleRequest();
 *
 * @copyright Genome Research Limited 2017
 */

const assert      = require('assert');
const send        = require('send');
const url         = require('url');
const LOGGER      = require('winston');
const DataMapper  = require('./mapper.js');
const DataAccess  = require('./auth.js');
const RangerModel = require('./model.js');
const HttpError   = require('./http/error.js');
const trailer     = require('./http/trailer.js');
const config      = require('../config.js');
const constants   = require('../constants.js');

const RANGER_PROXY_URL_KEY = 'rangerProxyURL';
const INVALID_INPUT        = 'InvalidInput';
const UNSUPPORTED_FORMAT   = 'UnsupportedFormat';

// Define file formats for send module

send.mime.default_type = 'application/octet-stream';
send.mime.define({
  'application/vnd.ga4gh.bam': ['bam'],
  'text/vnd.ga4gh.sam': ['sam'],
  'application/vnd.ga4gh.cram': ['cram'],
  'text/vnd.ga4gh.vcf': ['vcf']
});

/** Application controller */
class RangerController {

  /**
   * Creates a RangerController object instance.
   * @param request  - HTTP request object, required
   * @param response - HTTP response object, required
   * @param db       - Mongo db handle, required
   */
  constructor(request, response, db) {

    let validate = (obj) => {
      return (obj && (typeof obj == 'object'));
    };

    assert(validate(request), 'HTTP request object is required');
    this.request = request;

    assert(validate(response), 'Server response object is required');
    this.response = response;

    assert(validate(db), 'DB handle object is required');
    this.db = db;

    /**
     * Will this controller need to set a trailer as part of response.
     * @type {Boolean}
     */
    this.sendTrailer = trailer.trailersRequested(request.headers);

    this._user = null;
  }

  contentType(format) {
    let primaryContentType = RangerModel.isTextualFormat(format) ? 'text' : 'application';
    return [primaryContentType, 'vnd.ga4gh.' + format.toLowerCase()].join('/');
  }

  errorResponse(response, code, m, reasonPhrase) {
    // Create error object anyway, it will validate the input
    let err = new HttpError(response, code, m, reasonPhrase);
    if (!response.finished) {
      err.setErrorResponse();
      response.end();
    } else {
      throw new Error('Attempt to set error response on a closed stream');
    }
  }

  _methodIsAllowed() {
    let m = this.request.method;
    return m == 'GET' || m == 'OPTIONS';
  }

  _trimURL(url) {
    // Trim trailing forward slashes
    return url.replace(/\/*$/, '');
  }

  _setCORSHeaders() {

    let options = config.provide();
    this.response.setHeader('Vary', 'Origin');

    let requestOrigin = this.request.headers.origin;
    if (requestOrigin) {
      requestOrigin = this._trimURL(requestOrigin);
      let allowedOrigin = '';
      if (options.get('anyorigin')) {
        allowedOrigin = '*';
      } else if (options.get('originlist') &&
                 options.get('originlist').indexOf(requestOrigin) >= 0) {
        allowedOrigin = requestOrigin;
      }

      if (allowedOrigin) {
        this.response.setHeader('Access-Control-Allow-Origin', allowedOrigin);
        this.response.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS');
        this.response.setHeader('Access-Control-Allow-Headers', 'TE,X-Remote-User,withcredentials');
        this.response.setHeader('Access-Control-Max-Age', '1800'); // 30 minutes
        if (!options.get('skipauth')) {
          this.response.setHeader('Access-Control-Allow-Credentials', 'true');
        }
      }
    }
  }

  _sendData(query) {
    assert(query);
    assert(query.files.length > 0);

    if (query.files.length > 1) {
      this.errorResponse(this.response, 400,
        'More than one file exists with this name. ' +
        'Add a \'directory\' query parameter, then try again.',
        INVALID_INPUT);
      return;
    }
    if (query.files[0].startsWith('irods://')) {
      // file is on irods, we can't access it so fall back to using
      // the pipeline (samtools etc.)
      this._setupPipeline(query);
      return;
    }

    this._prepareTrailers();
    send(this.request, query.files[0]).pipe(this.response);
  }

  _prepareTrailers() {
    if (this.sendTrailer) {
      trailer.declare(this.response);
    }
  }

  _setupPipeline(query) {
    assert(query);

    this._prepareTrailers();
    this.response.setHeader('Content-Type', this.contentType(query.format));

    let endResponse = (truncated) => {
      if (this.sendTrailer) {
        trailer.setDataTruncation(this.response, truncated);
      }
      this.response.end();
    };

    try {
      let model = new RangerModel();
      model.process(query, this.response, endResponse);
    } catch (ex) {
      if (!ex.toString().startsWith('Inconsistent format')) {
        throw ex;
      }
      this.errorResponse(this.response, 500, ex);
    }
  }

  _returnReference(accession, reference) {
    let res = {
      accession: accession,
      reference: reference
    };
    let jsonRes = JSON.stringify(res);
    this.response.removeHeader('Transfer-Encoding');
    this.response.writeHead(
      200, 'OK', {
        'Content-Length': jsonRes.length,
        'Content-Type':   'application/json'});
    this.response.write(jsonRes);
    this.response.end();
  }

  _setDefaultListeners( dm ) {
    dm.once('error', ( err ) => {
      dm.removeAllListeners();
      this.errorResponse(this.response, 500, err);
      throw new Error(err);
    });
    dm.once('nodata', ( message ) => {
      dm.removeAllListeners();
      this.errorResponse(this.response, 404, message);
    });
  }

  _getReference( accession,  query ) {
    let options = config.provide();
    let host    = options.get('hostname');
    assert( accession, 'accession is required' );
    assert( host,      'host is required' );
    assert( query,     'query is required' );

    let dm = new DataMapper(this.db);
    this._setDefaultListeners(dm);
    dm.once('data', ( data ) => {
      dm.removeAllListeners();
      let reference = data[0].reference; // TODO index from query?
      LOGGER.debug(`Found reference ${reference} for acc ${query.accession}`);
      this._returnReference( accession, reference );
    });
    query.accession = accession;
    dm.getFileInfo(query, host);
  }

  _getData(query, authtype, cback) {
    let user = this._user;
    let options = config.provide();
    let host = options.get('hostname');
    assert((user && host && query));

    var dm = new DataMapper(this.db);
    this._setDefaultListeners(dm);
    dm.once('data', (data) => {
      query.files = data.map( (d) => {return d.file;} );
      query.reference = data[0].reference;
      LOGGER.debug(`Using ref file ${query.reference}`);
      dm.removeAllListeners();
      if (options.get('skipauth')) {
        cback(query);
      } else {
        let da = new DataAccess(authtype);
        // TODO unknown/undefined authtype should be caught in handleRequest(),
        // is an else clause needed?
        da.once('authorised', (username) => {
          da.removeAllListeners();
          LOGGER.debug(`User ${username} is given access`);
          cback(query);
        });
        da.once('failed', (username, message) => {
          da.removeAllListeners();
          LOGGER.info(message);
          this.errorResponse(this.response, 403,
            `Authorization failed for user '${username}'`);
        });
        let identifier;
        if (authtype === constants.AUTH_TYPE_USER) {
          identifier = user.username;
        } else if (authtype === constants.AUTH_TYPE_TOKEN) {
          identifier = user.token;
        }
        da.authorise(identifier, data.map( (d) => {return d.accessGroup;} ));
      }
    });
    dm.getFileInfo(query, host);
  }

  _checkFormat(format, requestFormat) {
    if ( !( format && requestFormat ) ) {
      throw new RangeError(`format and requestFormat must be defined`);
    }
    if (!RangerModel.supportsFormat(format)) {
      this.errorResponse(this.response, 400,
        `Format '${requestFormat}' is not supported, supported formats: ` +
        RangerModel.supportedFormats().join(', '),
        UNSUPPORTED_FORMAT
      );
      return false;
    }
    return true;
  }

  _ga4ghToDirectUrl(id, query) {
    LOGGER.debug('Process direct');
    LOGGER.debug(query);
    var _getInt = (str) => {
      /* Suppress unhelpful advice not to use Number as a constructor */
      /* jshint -W053 */
      let n = new Number(str);
      /* jshint +W053 */
      n = n.valueOf();
      if (!Number.isInteger(n)) {
        throw new RangeError(`'${str}' is not an integer`);
      }
      if (n < 0) {
        throw new RangeError(`'${str}' is not an unsigned integer`);
      }
      return n;
    };

    let url = `/sample?accession=${id}&format=${query.format}`;

    // Note that the returned 'internal' url will not feature any
    // unknown query parameters

    DataMapper.getFilterNames().forEach(function(parameter) {
      if ( parameter in query ) {
        url += '&' + parameter + '=';
        url += (query[parameter] !== undefined ? query[parameter] : '');
      }
    });

    if (query.referenceName) {
      let refSeparator = encodeURIComponent(':');
      let region = query.referenceName;
      assert(region.match(/^[-\w\.]+$/),
        'Invalid character in reference name ' + query.referenceName);
      /*
       * GA4GH ranges are semi-open starting from zero.
       * Samtools ranges are closed starting from one.
       */
      let start = query.start || '0';
      start = _getInt(start);
      let end = (typeof query.end === 'undefined') ? null : _getInt(query.end);
      if (end !== null && end <= start) {
        throw new RangeError('Range end should be bigger than start');
      }
      start += 1;

      if (end) {
        region += (refSeparator + start + '-' + end);
      } else {
        if (start !== 1) {
          region += refSeparator + start;
        }
      }
      url += ('&region=' + region);
    } else {
      if ((typeof query.start !== 'undefined') || (typeof query.end !== 'undefined')) {
        throw new Error(
          "'referenceName' attribute requered if 'start' or 'end' attribute is given");
      }
    }

    return url;
  }

  _serverUrlFromRequest() {
    let urlObj = url.parse(this.request.url, false);
    ['query', 'search', 'path', 'pathname', 'href']
      .forEach((key) => { urlObj[key] = undefined;});
    urlObj.slashes  = true;
    // Server URL is impossible to get from anywhere else
    urlObj.host     = this.request.headers.host;
    urlObj.protocol = config.provide().get('protocol');
    return url.format(urlObj);
  }

  _redirect(id, query) {
    let url = '';
    LOGGER.debug('Process redirect');
    LOGGER.debug(query);
    try {
      url = this._ga4ghToDirectUrl(id, query);
      LOGGER.debug('Forward url ' + url);
    } catch (err) {
      this.errorResponse(this.response, 400, err.message, INVALID_INPUT);
    }
    if (url) {
      let serverUrl = query[RANGER_PROXY_URL_KEY] ?
        query[RANGER_PROXY_URL_KEY] :
        this._serverUrlFromRequest();
      let reply = {};
      reply.urls   = [{url: serverUrl + url}];
      reply.format = query.format;
      let redirect = JSON.stringify(reply);
      LOGGER.debug('JSON produced ' + redirect);
      this.response.removeHeader('Transfer-Encoding');
      this.response.writeHead(
        200, 'OK, see redirection instructions in the body of the message', {
          'Content-Length': redirect.length,
          'Content-Type':   'application/json'});
      this.response.write(redirect);
      this.response.end();
    }
  }

  _fileClback(q, authtype) {
    if (!q.name) {
      this.errorResponse(this.response, 400,
        'Invalid request: file name should be given',
        INVALID_INPUT
      );
    } else {
      this._getData(q, authtype, this._sendData.bind(this));
    }
  }


  _sampleClback(q, authtype) {
    if (!q.accession) {
      this.errorResponse(this.response, 400,
        'Invalid request: sample accession number should be given',
        INVALID_INPUT
      );
    } else {
      this._getData(q, authtype, this._setupPipeline.bind(this));
    }
  }

  _redirectClback(q, authtype, p, id) {
    if (!id) {
      this.errorResponse(this.response, 400,
        `Invalid request ${p}: sample accession number should be given`,
        INVALID_INPUT
      );
    } else {
      this._redirect(id, q);
    }
  }

  _referenceClback(q, authtype, p, id) {
    if (!id) {
      this.errorResponse(this.response, 400,
        `Invalid request ${p}: sample accession number should be given`,
        INVALID_INPUT
      );
    } else {
      this._getReference(id, q);
    }
  }

  _detectProxy(authtype) {
    let knownProxies = config.provide().get('proxylist');
    if (knownProxies) {
      let fhost = this.request.headers['x-forwarded-host'];
      if (!fhost) {
        throw new Error('Bypassing proxy server is not allowed');
      }
      let fproto = config.provide().get('protocol');
      fhost = fproto + '//' + fhost;
      fhost = this._trimURL(fhost);
      if (Object.keys(knownProxies).indexOf(fhost) < 0) {
        throw new Error(`Unknown proxy ${fhost}`);
      }
      // If the ga4gh url was requested using an authorisation type, we must
      // provide a matching authtype redirect in the returned url.
      return knownProxies[fhost] + (authtype ? ('/' + authtype) : '');
    }
    return null;
  }

  /**
   * <p>Interprets a query url and its parameters.
   * Maps the query to data sources and authorises the client for access.
   * Handles request and sends a reply to the client.</p>
   *
   * <p>Only GET and OPTIONS requests are supported.</p>
   *
   * <p>User authentication should be performed elsewhere. It is expected that the
   * 'x-remote-user' header is set in the request.</p>
   *
   * <p>User authorisation is per data resource, see 'auth' module for details.</p>
   *
   * <p>The following URLs are accepted:</p>
   * <pre>
   *   URL           Paramaters
   *   /file         name      - required
   *                 directory - optional
   *                 format    - required data format; 'SAM', 'CRAM', 'BAM' and
   *                             'VCF' are recognised; if the format is not one of
   *                             these four, error response is returned; if the
   *                             format is not given, 'BAM' is used
   *                 region    - optional, in the format 'chr1:3000-3020', multiple
   *                             key-value pairs are allowed
   *
   *   /sample       accession - sample accession number, required
   *                 region    - as for /file, but a single key-value pair only
   *                 format    - as for /file
   *
   *   /sample/XXXX/reference
   *     where XXXX is a sample accession number (the same as accession
   *     parameter in the /sample url)
   *
   *   /ga4gh/v.0.1/get/sample/XXXX
   *     where XXXX is a sample accession number (the same as accession
   *     parameter in the /sample url)
   *                 format    - as above
   *                 referenceName - chromosome name, optional
   *                 start     - The start position of the range on the reference,
   *                             0-based inclusive; optional, if specified, referenceName
   *                             must also be specified
   *                 end       - The end position of the range on the reference,
   *                             0-based exclusive; optional, if specified, referenceName
   *                             must also be specified
   * </pre>
   *
   * <p>The URL starting with /ga4gh/v.0.1 outputs redirection JSON response according to
   * to the version of the GA4GH specification given in the URL. The rest of URLs stream the
   * data directly.</p>
   *
   * <p>Redirection JSON response contains a single absolute URL. In production this server is
   * likely to be running behind a reverse proxy. In this situation the JSON redirection respose
   * will contain the URL on the proxy server.</p>
   *
   * <p>The controller uses X-Forwarded-Host request header to determine the URL
   * of the proxy server.</p>
   */
  handleRequest() {

    let options = config.provide();

    if (!this._methodIsAllowed()) {
      return this.errorResponse(
        this.response, 405, this.request.method + ' method is not allowed');
    }

    let urlObj = url.parse(this.request.url, true);
    let path = urlObj.pathname || '';
    let q = urlObj.query;

    // RegExp below is equivalent to /^\/(authuser|authtoken)/
    let authmatch = path.match(new RegExp('^\/(' + constants.AUTH_TYPE_USER + '|' + constants.AUTH_TYPE_TOKEN + ')'));
    let authtype = authmatch ? authmatch[1] : null;

    // authtype has been found, and the rest of the controller does not expect to see `/${authtype}` in url, so strip
    if (authmatch) {
      path = path.replace(authmatch[0], '');
    }

    let proxyURL;
    try {
      proxyURL = this._detectProxy(authtype);
    } catch (err) {
      this.errorResponse(this.response, 403, err.message);
    }
    if (typeof proxyURL === 'undefined') {
      return;
    }

    this._setCORSHeaders();
    if (this.request.method == 'OPTIONS') {
      this.response.end();
      return;
    }

    this._user = {};
    if (!options.get('skipauth')) {
      if (authtype === constants.AUTH_TYPE_USER) {
        this._user.username = this.request.headers['x-remote-user'] || null;
        if (!this._user.username) {
          return this.errorResponse(this.response, 401, 'Proxy authentication required');
        }
      } else if (authtype === constants.AUTH_TYPE_TOKEN) {
        this._user.token = this.request.headers[constants.TOKEN_BEARER_KEY_NAME] || null;
        if (!this._user.token) {
          return this.errorResponse(this.response, 401, `Bearer token required; visit ${options.get('authurl')}`);
        }
      } else if (!authtype) {
        return this.errorResponse(this.response, 401,
          'authorisation type must be specified in url');
      } else {
        return this.errorResponse(this.response, 401, 'Unrecognised authtype');
      }
    }

    let duplicateAttr = (query) => {
      return Object.keys(query).find(
        (element) => { return typeof query[element] === 'object'; }
      );
    };
    // We disallow query strings containing multiple attribute-value pairs
    // for the same attribute. In most places our code expects attribute
    // values to be strings and breaks if the value is an array.
    // We might opt for more granular consideration later.
    let dup = duplicateAttr(q);
    if (dup) {
      return this.errorResponse(this.response, 400,
        `Invalid request: multiple values for attribute '${dup}'`,
        INVALID_INPUT
      );
    }

    q[RANGER_PROXY_URL_KEY] = proxyURL;

    let requestFormat;
    if (!q.format) {
      q.format = requestFormat = RangerModel.defaultFormat();
    } else {
      requestFormat = q.format;
      q.format = q.format.toUpperCase();
    }
    if (!this._checkFormat(q.format, requestFormat)) {
      return;
    }
    if (q.format === 'VCF' && options.get('multiref')) {
      return this.errorResponse(this.response, 400,
        'Invalid request: ' +
        'cannot produce VCF files while multiref set on server',
        INVALID_INPUT
      );
    }

    var matchPath = (p) => {
      let map = new Map();
      map.set(new RegExp(/^\/file$/),   'fileClback');
      map.set(new RegExp(/^\/sample$/), 'sampleClback');
      map.set(new RegExp(/^\/ga4gh\/v\.0\.1\/get\/sample\/([\.\w]+)$/), 'redirectClback');
      map.set(new RegExp(/^\/sample\/([\.\w]+)\/reference$/), 'referenceClback');

      let it = map.entries();
      let el = it.next();
      let reResult = null;
      let callback = null;
      while (!reResult && el && el.done === false) {
        let re = el.value[0];
        reResult = re.exec(p);
        callback = el.value[1];
        el = it.next();
      }
      return !reResult ? false : {path: reResult[0], captured: reResult[1], cback: callback};
    };

    let match = matchPath(path);
    if (!match) {
      this.errorResponse(this.response, 404, 'URL not found : ' + path);
    } else {
      let cback = '_' + match.cback;
      assert(typeof this[cback] === 'function', 'Callback function is not defined');
      this[cback](q, authtype, match.path, match.captured);
    }
  }
}

module.exports = RangerController;
