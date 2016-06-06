"use strict";

/**
 * controller module.
 * @module lib/controller
 *
 * @description Controller module for the streaming HTTP server.
 *
 * @example <caption>Example usage of the controller.</caption>
 *   const RangerController = require('../lib/controller.js');
 *   // Create a new object, no authorisation.
 *   let c = new RangerController(request, response, db);
 *   // Create a new object with authorisation.
 *   let c1 = new RangerController(request, response, db, dir, true);
 *   // handle request
 *   c.handleRequest(host);
 *
 * @author Marina Gourtovaia
 * @copyright Genome Research Limited 2016
 */

const assert      = require('assert');
const fs          = require('fs');
const http        = require('http');
const os          = require('os');
const url         = require('url');
const DataAccess  = require('./auth.js');
const DataMapper  = require('./mapper.js');
const RangerModel = require('./model.js');
const HttpError   = require('./http/error.js');
const trailer     = require('./http/trailer.js');

/** Application controller */
class RangerController {

  /**
   * Creates a RangerController object instance.
   * @param request  - HTTP request object, required
   * @param response - HTTP response object, required
   * @param db       - Mongo db handle, required
   * @param tmpDir   - directory path for temporary data, optional,
   *                   default to an OS tmp directory
   * @param skipAuth - boolen flag, optional, if true instructs to
   *                   skip authorisation
   */
  constructor(request, response, db, tmpDir, skipAuth) {

    var validate = (obj) => {
      return (obj && (typeof obj == 'object'));
    };

    assert(validate(request), 'HTTP request object is required');
    this.request = request;

    assert(validate(response), 'Server response object is required');
    assert((response instanceof http.ServerResponse),
      'Server response object is required');
    this.response = response;

    assert(validate(db), 'DB handle object is required');
    this.db = db;

    this.tmpDir = (tmpDir || os.tmpdir());
    assert(fs.existsSync(this.tmpDir),
      `Temp data directory '${this.tmpDir}' does not exist`);
    this.skipAuth = !!skipAuth;
  }

  contentType(format) {
    let primaryContentType = RangerModel.isTextualFormat(format) ? 'text' : 'application';
    return [primaryContentType, format].join('/');
  }

  errorResponse(response, code, m) {
    // Create eror object anyway, it will validate the input
    let err = new HttpError(response, code, m, true);
    if (!response.finished) {
      err.setErrorResponse();
      response.end();
    } else {
      throw new Error('Attempt to set error response on a closed stream');
    }
  }

  _setupPipeline(query) {
    assert(query);

    trailer.declare(this.response);
    this.response.setHeader("Content-Type", this.contentType(query.format));

    var endResponse = (truncated) => {
      trailer.setDataTruncation(this.response, truncated);
      this.response.end();
    };

    try {
      let model = new RangerModel(this.tmpDir);
      model.process(query, this.response, endResponse);
    } catch (ex) {
      if (!ex.toString().startsWith('Inconsistent format')) {
        throw ex;
      }
      this.errorResponse(this.response, 500, ex);
    }
  }

  _getData(user, host, query) {
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
      dm.removeAllListeners();
      if (this.skipAuth) {
        this._setupPipeline(query);
      } else {
        var da = new DataAccess(this.db);
        da.once('authorised', (username) => {
          da.removeAllListeners();
          console.log(`User ${username} is given access`);
          this._setupPipeline(query);
        });
        da.once('failed', (username, message) => {
          da.removeAllListeners();
          this.errorResponse(this.response, 403,
            `Authorisation failed for user '${username}': ${message}`);
        });
        da.authorise(user.username, data.map( (d) => {return d.accessGroup;} ));
      }
    });
    dm.getFileInfo(query, host);
  }

  _checkFormat(format) {
    if (!RangerModel.supportsFormat(format)) {
      this.errorResponse(this.response, 409,
        `Format '${format}' is not supported, supported formats: ` +
        RangerModel.supportedFormats().join(', '));
      return false;
    }
    return true;
  }

  _ga4ghToDirectUrl(id, query) {

    var _getInt = (str) => {
      /* Suppress unhelpful advice not to use Number as a constructor */
      /* jshint -W053 */
      let n = new Number(str);
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
      let start = query.start || '0';
      start = _getInt(start) + 1;
      if (query.end) {
        let end = _getInt(query.end) + 1;
        if (end < start) {
          throw new RangeError('Range end should be bigger that start');
        }
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

  _redirect(id, query) {
    let url = '';
    try {
      url = this._ga4ghToDirectUrl(id, query);
    } catch (err) {
      this.errorResponse(this.response, 422, err.message);
    }
    if (url) {
      let reply = {};
      reply.urls   = [url];
      reply.format = query.format.toUpperCase();
      let redirect = JSON.stringify(reply);
      this.response.removeHeader('Transfer-Encoding');
      this.response.writeHead(200, 'OK, see redirection instructions in the body of the message', {
        'Content-Length': redirect.length,
        'Content-Type':   'application/json'});
      this.response.write(redirect);
      this.response.end();
    }
  }

  /**
   * Interprets a query url and its parameters.
   * Maps the query to data sources and authorises the client for access.
   * Handles request and sends a reply to the client.
   *
   * User authentication should be performed elsewhere. It is expected that the
   * 'x-remote-user' header is set in the request.
   *
   * User authorisation is per data resource, see 'auth' module for details.
   *
   * @param host - the name of the host to retrieve the data from.
   *               If the data are not available on the given host, the
   *               default host will be used.
   *
   * The following URLs are accepted:
   *   URL           Paramaters
   *   /file         name      - required
   *                 directory - optional
   *                 format    - required data format; 'sam', 'cram' and 'bam'
   *                             are recognised; if the format is not one of
   *                             these three, error response is returned; if the
   *                             format is not given, 'bam' is used
   *                 region    - optional, in the format 'chr1:3000-3020', multiple
   *                             key-value pairs are allowed
   *   /sample       accession - sample accession number, required
   *                 region    - as for /file, but a single key-value pair only
   *                 format    - as for /file
   *   /api/ga4gh/v.0.1/get/sample/XXXX, where XXXX is a sample accession number
   *                 format    - as above
   *                 referenceName - chromosome name, optional
   *                 start     - The start position of the range on the reference,
   *                             0-based inclusive; optional, if specified, referenceName
   *                             must also be specified
   *                 end       - The end position of the range on the reference,
   *                             0-based exclusive; optional, if specified, referenceName
   *                             must also be specified
   */
  handleRequest(host) {
    assert(host, 'The data host name is required');

    let user = {};
    user.username   = this.request.headers['x-remote-user'] || null;
    if (!user.username && !this.skipAuth) {
      return this.errorResponse(this.response, 401, 'Proxy authentication required');
    }

    let urlObj = url.parse(this.request.url, true);
    let path = urlObj.pathname;
    path = path ? path : '';
    let q = urlObj.query;

    if (!q.format) {
      q.format = RangerModel.defaultFormat();
    }
    if (!this._checkFormat(q.format)) {
      return;
    }

    let directClback = (p) =>  {
      if (p === '/file' && !q.name) {
        this.errorResponse(this.response, 422,
          'Invalid request: file name should be given');
      } else if (p === '/sample' && !q.accession) {
        this.errorResponse(this.response, 422,
          'Invalid request: sample accession number should be given');
      } else {
        this._getData(user, host, q);
      }
    };

    let redirectClback = (p, id) =>  {
      if (!id) {
        this.errorResponse(this.response, 422,
          `Invalid request ${p}: sample accession number should be given`);
      } else {
        this._redirect(id, q);
      }
    };

    var matchPath = (p) => {
      let map = new Map();
      map.set(new RegExp(/^\/file$/),   directClback);
      map.set(new RegExp(/^\/sample$/), directClback);
      map.set(new RegExp(/^\/api\/ga4gh\/v\.0\.1\/get\/sample\/([\.\w]+)$/), redirectClback);

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
      let cback = match.cback;
      cback(match.path, match.captured);
    }
  }
}

module.exports = RangerController;