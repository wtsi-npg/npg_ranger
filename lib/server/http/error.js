"use strict";

/**
 * error module.
 * @module server/http/error
 *
 * @description GA4GH API-compliant HTTP errors.
 *
 * @example <caption>Example usage of the error module.</caption>
 *   const ServerHttpError = require('../lib/server/http/error.js');
 *   const HttpError       = ServerHttpError.HttpError;
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

const HTTP_ERROR_TYPES = {
  400: 'BadRequest',
  401: "InvalidAuthentication",
  403: "PermissionDenied",
  404: "NotFound",
  405: "MethodNotAllowed",
  500: "InternalError",
  503: "ServiceUnavailable",
};

const INVALID_INPUT      = 'InvalidInput';
const UNSUPPORTED_FORMAT = 'UnsupportedFormat';

/** A helper class to set HTTP error response */
class HttpError {
  /**
   * Creates an HttpError type object. Validates input. Replaces error code for
   * 500 if an unknown error code is used..
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

    let codeString = '' + code;

    if ( HTTP_ERROR_TYPES.hasOwnProperty(codeString) ) {
      this.code      = code;
      this.errorType = ( reasonPhrase ) ? reasonPhrase
                                        : HTTP_ERROR_TYPES[codeString];
      this.message   = ( message || 'Unknown error' );
      LOGGER.info(`Server error ${code}: ${this.errorType} | ${this.message}.`);
    } else {
      LOGGER.error(`Failed to process unknown error code: '${code}', message: '${message}'`);
      this.code      = 500;
      this.errorType = HTTP_ERROR_TYPES[500];
      this.message   = 'Internal server error';
    }
  }

  /**
   * Sets the content type to 'application/json', response code as passed to
   * the object's constructor and the body of the response to a json string
   * according to the GA4GH specification. Chunked encoding is disabled.
   *
   */
  setErrorResponse() {
    let jerr = JSON.stringify({
      htsget: {
        error:   this.errorType,
        message: this.message
      }
    });
    this.response.removeHeader('Transfer-Encoding');
    this.response.writeHead(this.code, this.message, {
      'Content-Length': jerr.length,
      'Content-Type':   'application/json'});
    this.response.write(jerr);
  }
}

module.exports = {
  HttpError,
  INVALID_INPUT,
  UNSUPPORTED_FORMAT
};
