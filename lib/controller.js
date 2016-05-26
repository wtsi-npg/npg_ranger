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
const os          = require('os');
const url         = require('url');
const DataAccess  = require('./auth.js');
const DataMapper  = require('./mapper.js');
const RangerModel = require('./model.js');
const HttpError   = require('./http/error.js');
const trailer     = require('./http/trailer.js');

const DEFAULT_FORMAT  = 'bam';

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

    assert(request, 'HTTP request object is required');
    this.request = request;

    assert(response, 'HTTP response object is required');
    this.response = response;

    assert(db, 'DB handle is required');
    this.db = db;

    this.tmpDir = (tmpDir || os.tmpdir());
    this.skipAuth = skipAuth;
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
    if (!query.format) {
      query.format = DEFAULT_FORMAT;
    }
    if (query.format === 'bam' || query.format === 'cram') {
      this.response.setHeader("Content-Type", 'application/octet-stream');
    }

    var endResponse = (truncated) => {
      trailer.setDataTruncation(this.response, truncated);
      this.response.end();
    };
    let model = new RangerModel(this.response,
                                query,
                                this.tmpDir);
    if (query.files.length === 1) {
      model.getFile(endResponse);
    } else {
      try {
        model.mergeFiles(endResponse);
      } catch (ex) {
        if (!ex.toString().startsWith('Inconsistent format')) {
          throw ex;
        }
        this.errorResponse(this.response, 500, ex);
      }
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
   *                             these three or not given, 'bam' is used
   *                 range     - optional, in the format 'chr13000:3020', multiple
   *                             key-value pairs are allowed
   *   /sample       accession - sample accession number, required
   *                 range     - as for /file, but a single key-value pair only
   *                 format    - as for /file
   */
  handleRequest(host) {
    assert(host);

    let user = {};
    user.username   = this.request.headers['x-remote-user'] || null;
    if (!user.username && !this.skipPath) {
      return this.errorResponse(this.response, 401, 'Proxy authentication required');
    }

    let url_obj = url.parse(this.request.url, true);
    let path = url_obj.pathname;
    path = path ? path : '';
    let q = url_obj.query;

    switch (path) {
      case '/file': {
        if (!q.name) {
          this.errorResponse(this.response, 422,
            'Invalid request: file name should be given');
        } else {
          this._getData(user, host, q);
        }
        break;
      }
      case '/sample': {
        if (!q.accession) {
          this.errorResponse(this.response, 422,
            'Invalid request: sample accession number should be given');
        } else {
          this._getData(user, host, q);
        }
        break;
      }
      default: {
        this.errorResponse(this.response, 422, 'URL not available: ' + path);
      }
    }
  }
}

module.exports = RangerController;