"use strict";

/**
 * error module.
 * @module server/http/error
 *
 * @description GA4GH API-compliant HTTP errors.
 *
 * @example <caption>Example usage of the error module.</caption>
 *   const HttpError = require('../lib/server/http/error.js');
 *   // Create a new object.
 *   let err = new HttpError(response, code, m);
 *   // Create a new object with tailored status phrase
 *   let err = new HttpError(response, code, m, myStatusPhrase);
 *   // Set response
 *   err.setErrorResponse();
 *
 * @copyright Genome Research Limited 2017
 */

const LOGGER  = require('winston');

const trailer = require('./trailer.js');

const HTTP_ERROR_TYPES = {
  400: 'BadRequest',
  401: "InvalidAuthentication",
  403: "PermissionDenied",
  404: "NotFound",
  405: "MethodNotAllowed",
  500: "InternalError",
  503: "ServiceUnavailable",
};

/** A helper class to set HTTP error response */
class HttpError {
  /**
   * Creates an HttpError type object. Validates input. Raises
   * an error if the error code in not one of valid for this
   * application error codes.
   * @param response     - HTTP response object
   * @param code         - HTTP error code
   * @param message      - error message, can be undefined
   * @param reasonPhrase - reason phrase for status if different from standard
   */
  constructor(response, code, message, reasonPhrase) {

    if (!response) {
      throw new ReferenceError('HTTP response object is required');
    }
    this.response = response;

    if (!code) {
      throw new ReferenceError('HTTP error code is required');
    }
    let errorType = ( reasonPhrase ) ? reasonPhrase
                                     : HTTP_ERROR_TYPES[code + '']; // code as string
    if ( errorType === undefined ) {
      throw new RangeError(`Invalid HTTP error code: ${code}`);
    }
    this.code      = code;
    this.errorType = errorType;

    let m = '' + (message || 'Unknown error');
    this.message   = m;
    LOGGER.info(`Server error ${code}: ${this.message}.`);
  }

  /**
   * If the headers have not been sent yet, sets the content type
   * to 'application/json', response code as passed to the object's
   * constructor and the body of the response to a json string
   * according to the GA4GH specification. Chunked encoding is
   * disabled.
   *
   * If the headers have been sent, the trailer is set to mark the
   * data truncated.
   */
  setErrorResponse() {
    if (!this.response.headersSent) {
      let jerr = JSON.stringify({
        error:   this.errorType,
        message: this.message
      });
      trailer.removeDeclaration(this.response);
      this.response.removeHeader('Transfer-Encoding');
      this.response.writeHead(this.code, this.message, {
        'Content-Length': jerr.length,
        'Content-Type':   'application/json'});
      this.response.write(jerr);
    } else {
      trailer.setDataTruncation(this.response, true);
    }
  }
}

module.exports = HttpError;
