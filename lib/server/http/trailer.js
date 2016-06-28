"use strict";

/**
 * trailer module.
 * @module lib/server/http/trailer
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

function validateResponse(response) {
  if (!response) {
    throw new ReferenceError('HTTP response object is required');
  }
}

/**
 * Declares the names of the trailers. If called after the headers
 * have been set, an error is raised.
 * @param response - HTTP response object
 */
exports.declare = (response) => {
  validateResponse(response);
  response.setHeader(TRAILER_HEADER_NAME, DATA_TRUNCATION_TRAILER);
};

/**
 * Removes trailer declaration
 * @param response - HTTP response object
 */
exports.removeDeclaration = (response) => {
  validateResponse(response);
  response.removeHeader(TRAILER_HEADER_NAME);
};

/**
 * Sets data truncation trailer to the given value.
 * @param response  - HTTP response object
 * @param truncated - boolean flag indicating whether the data
 *                    are truncated
 */
exports.setDataTruncation = (response, truncated) =>  {
  validateResponse(response);
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
