"use strict";

/**
 * trailer module.
 * @module server/http/trailer
 *
 * @description Helper module for setting HTTP trailer headers.
 *
 * @example <caption>Example usage of the trailer module.</caption>
 *   const trailer = require('../lib/server/http/trailer.js');
 *   // Declare trailers
 *   trailer.declare(response);
 *   // Mark the data as truncated
 *   trailer.setDataTruncation(response);    // Data truncated
 *   trailer.setDataTruncation(response, 0); // Data truncated
 *   trailer.setDataTruncation(response, 1); // Data not truncated
 *
 * @author Marina Gourtovaia
 * @copyright Genome Research Limited 2016
 */

const TRAILER_HEADER_NAME     = 'Trailer';
const DATA_TRUNCATION_TRAILER = 'data-truncated';

function validate(object, type) {
  type = type || 'response';
  if (!object) {
    throw new ReferenceError(`HTTP ${type} object is required`);
  }
}

/**
 * Declares the names of the trailers. If called after the headers
 * have been set, an error is raised.
 * @param response - HTTP response object
 */
exports.declare = (response) => {
  validate(response);
  response.setHeader(TRAILER_HEADER_NAME, DATA_TRUNCATION_TRAILER);
};

/**
 * Removes trailer declaration
 * @param response - HTTP response object
 */
exports.removeDeclaration = (response) => {
  validate(response);
  response.removeHeader(TRAILER_HEADER_NAME);
};

/**
 * Sets data truncation trailer to the given value.
 * @param response  - HTTP response object
 * @param truncated - boolean flag indicating whether the data
 *                    are truncated
 */
exports.setDataTruncation = (response, truncated) =>  {
  validate(response);
  if ( typeof (truncated) !== 'boolean' ) {
    throw new ReferenceError('boolean flag indicating data truncation is required');
  }
  let trailers = response.getHeader(TRAILER_HEADER_NAME);
  // Check that this trailer has been declared
  if (trailers && trailers.includes(DATA_TRUNCATION_TRAILER)) {
    var header = {};
    header[DATA_TRUNCATION_TRAILER] = truncated.toString();
    response.addTrailers(header);
  } else {
    throw new Error('Cannot set data truncation trailer because it has not been declared');
  }
};

function trailersRequested(request) {
  validate(request, 'request');
  return request.headers &&
         request.headers.te &&
         request.headers.te.includes('trailers');
}

/**
 * Returns true if the request states that the client can accept
 * trailers.
 * @param request - HTTP request object
 */
exports.trailersRequested = trailersRequested;

/**
 * Returns true if the request stated that the client can accept
 * trailers and response contains a trailer indicating that
 * that the response data is not fit for purpose (truncated).
 * @param request  - HTTP request object
 * @param response - HTTP response object
 */
exports.isDataTruncated = (request, response) => {
  validate(request, 'request');
  validate(response);
  return trailersRequested(request) &&
         response.trailers &&
         response.trailers[DATA_TRUNCATION_TRAILER] &&
         response.trailers[DATA_TRUNCATION_TRAILER] === true.toString();
};

