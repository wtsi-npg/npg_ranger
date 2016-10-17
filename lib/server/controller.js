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
 * @author Marina Gourtovaia
 * @copyright Genome Research Limited 2016
 */

const assert      = require('assert');
const url         = require('url');
const LOGGER      = require('winston');
const DataAccess  = require('./auth.js');
const DataMapper  = require('./mapper.js');
const RangerModel = require('./model.js');
const HttpError   = require('./http/error.js');
const trailer     = require('./http/trailer.js');
const config      = require('../config.js');

const RANGER_PROXY_URL_KEY = 'rangerProxyURL';

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

  errorResponse(response, code, m) {
    // Create error object anyway, it will validate the input
    let err = new HttpError(response, code, m, true);
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
        this.response.setHeader('Access-Control-Allow-Headers', 'TE,X-Remote-User');
        this.response.setHeader('Access-Control-Max-Age', '1800'); // 30 minutes
        if (!options.get('skipauth')) {
          this.response.setHeader('Access-Control-Allow-Credentials', 'true');
        }
      }
    }
  }

  _setupPipeline(query) {
    assert(query);

    if (this.sendTrailer) {
      trailer.declare(this.response);
    }
    this.response.setHeader('Content-Type', this.contentType(query.format));

    var endResponse = (truncated) => {
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

  _getData(query) {

    let user = this._user;
    let options = config.provide();
    let host = options.get('hostname');
    assert((user && host && query));

    var dm = new DataMapper(this.db);
    dm.once('error', (err) => {
      dm.removeAllListeners();
      this.errorResponse(this.response, 500, err);
    });
    dm.once('nodata', (message) => {
      dm.removeAllListeners();
      this.errorResponse(this.response, 404, message);
    });
    dm.once('data', (data) => {
      query.files = data.map( (d) => {return d.file;} );
      query.reference = data[0].reference;
      LOGGER.debug(`Using ref file ${query.reference}`);
      dm.removeAllListeners();
      if (options.get('skipauth')) {
        this._setupPipeline(query);
      } else {
        var da = new DataAccess(this.db);
        da.once('authorised', (username) => {
          da.removeAllListeners();
          LOGGER.debug(`User ${username} is given access`);
          this._setupPipeline(query);
        });
        da.once('failed', (username, message) => {
          da.removeAllListeners();
          LOGGER.info(message);
          this.errorResponse(this.response, 403,
            `Authorization failed for user '${username}'`);
        });
        da.authorise(user.username, data.map( (d) => {return d.accessGroup;} ));
      }
    });
    dm.getFileInfo(query, host);
  }

  _checkFormat(format, requestFormat) {
    if ( !( format && requestFormat ) ) {
      throw new RangeError(`format and requestFormat must be defined`);
    }
    if (!RangerModel.supportsFormat(format)) {
      this.errorResponse(this.response, 409,
        `Format '${requestFormat}' is not supported, supported formats: ` +
        RangerModel.supportedFormats().join(', '));
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
      this.errorResponse(this.response, 422, err.message);
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

  _directClback(q, p) {
    if (p === '/file' && !q.name) {
      this.errorResponse(this.response, 422,
        'Invalid request: file name should be given');
    } else if (p === '/sample' && !q.accession) {
      this.errorResponse(this.response, 422,
        'Invalid request: sample accession number should be given');
    } else {
      this._getData(q);
    }
  }

  _redirectClback(q, p, id) {
    if (!id) {
      this.errorResponse(this.response, 422,
        `Invalid request ${p}: sample accession number should be given`);
    } else {
      this._redirect(id, q);
    }
  }

  _detectProxy() {
    let knownProxies = config.provide().get('proxylist');
    if (knownProxies) {
      let fhost = this.request.headers['x-forwarded-host'];
      if (!fhost) {
        throw new Error('Bypassing proxy server is not allowed');
      }
      let fproto = this.request.headers['x-forwarded-proto'] ||
                   config.provide().get('protocol');
      fhost = fproto + '//' + fhost;
      fhost = this._trimURL(fhost);
      if (Object.keys(knownProxies).indexOf(fhost) < 0) {
        throw new Error(`Unknown proxy ${fhost}`);
      }
      return knownProxies[fhost];
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
   * <p>The controller uses X-Forwarded-Host request header to determine the URL of the proxy server.
   * The protocol is inferred from X-Forwarded-Proto request header and defaults to http. If the
   * X-Forwarded-Host-Suffix request header is available, its value will be appended to the
   * server url (no additional separators are used).</p>
   */
  handleRequest() {

    if (!this._methodIsAllowed()) {
      return this.errorResponse(
        this.response, 405, this.request.method + ' request is not allowed');
    }

    let proxyURL;
    try {
      proxyURL = this._detectProxy();
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

    let options = config.provide();
    this._user = {};
    this._user.username   = this.request.headers['x-remote-user'] || null;
    if (!this._user.username && !options.get('skipauth')) {
      return this.errorResponse(this.response, 401, 'Proxy authentication required');
    }

    let urlObj = url.parse(this.request.url, true);
    let path = urlObj.pathname;
    path = path ? path : '';
    let q = urlObj.query;
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
      return this.errorResponse(this.response, 422,
        'Invalid request: ' +
        'cannot produce VCF files while multiref set on server');
    }

    var matchPath = (p) => {
      let map = new Map();
      map.set(new RegExp(/^\/file$/),   'directClback');
      map.set(new RegExp(/^\/sample$/), 'directClback');
      map.set(new RegExp(/^\/ga4gh\/v\.0\.1\/get\/sample\/([\.\w]+)$/), 'redirectClback');

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
      this[cback](q, match.path, match.captured);
    }
  }
}

module.exports = RangerController;
