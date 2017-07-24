"use strict";

const assert = require('assert');

/**
 * NodeJS implementation of assert
 * @external assert
 * @see  {@link https://nodejs.org/dist/latest-v6.x/docs/api/assert.html|assert}
 */

/**
 * @module client/uriUtils
 *
 * @description
 * <p>Module provinding some basic utils to deal with URIs.</p>
 *
 * @requires {@link external:assert|assert}
 *
 * @copyright Genome Research Limited 2017
 */

/**
 * Builds a &lt;Buffer&gt; object from a data URI using the encoding provided as
 * part of the data URI. Decoding is delegated to the &lt;Buffer&gt; implementation.
 *
 * @example
 *
 * let encodedData = new Buffer('some data', 'ascii').toString('base64');
 * let dataURI = 'data:text/plain;charset=utf-8;base64,' + encodedData;
 * let buffer = uriUtils.procDataURI(dataURI);
 * process.stdout.write(buffer.toString() + '\n'); // some data
 *
 * @param  {string} uri data URI
 * @return {Buffer}     Buffer object created with defined encoding
 */
var procDataURI = ( uri ) => {
  assert(uri, 'uri is required');
  let regex = /^data:[\S]+\/[\S]+;([\S]+),(.*)/;

  let matches = uri.match(regex);
  let buffer;

  if ( matches ) {
    let encoding = matches[1];
    let data     = matches[2];

    buffer = new Buffer(data, encoding);
  } else {
    throw new Error('Unable to decode uri: ' + uri);
  }

  return buffer;
};

module.exports = {
  procDataURI: procDataURI
};
