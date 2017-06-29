"use strict";

const assert = require('assert');
const util   = require('util');

const TOKEN_REGEX          = /^\s*Bearer\s+(\S+)\s*$/i;
const TOKEN_BEARER_FORMAT  = `Bearer %s`;
const TOKEN_BEARER_MAX_LEN = 1000;
/**
 * Used for testing for invalid characters in the Authorization header value
 * @type {RegExp}
 */
const TOKEN_VALUE_REGEX    = /^[a-zA-Z0-9\+\/\-\_ =]+$/;


/**
 * Used by client to format a string token to be used as an Authorization header.
 *
 * @param  {String} token the string with the original token
 * @return {String}       a string with the token formatted for use as Auth
 *                        header. It will match 'Bearer <token>'
 * @throws Error when an empty (zero length or blank) token is passed as parameter
 */
let formatTokenForHeader = token => {
  assert(typeof token === 'string', 'String parameter is required');
  token = token.trim();
  assert(token.length > 0, 'token cannot be blank');
  return util.format(TOKEN_BEARER_FORMAT, token);
};

/**
 * Used by server to parse a token from the 'Authorization' header. It is expected
 * the string will match 'Bearer <token>' pattern.
 * @param  {String} string value of the original header
 * @return {String}        token parsed from header value
 * @throws Error when the string parameter is contains unexpected characters
 *         alphanumeric and </+-_>
 * @throws Error when the string parameter is too long
 * @throws Error when the string does not match the expected pattern 'Bearer <token>'
 */
let parseToken = string => {
  assert(typeof string === 'string', 'String parameter is required');
  if ( string.length > 1000 ) {
    throw new Error(
      `String provided for parsing token is too long (>${TOKEN_BEARER_MAX_LEN})`
    );
  }
  if ( !TOKEN_VALUE_REGEX.test(string) ) {
    throw new Error('Token contains unexpected characters');
  }
  let parseResult = TOKEN_REGEX.exec(string);
  if ( parseResult === null ) {
    throw new Error('Unexpected format in authorization string.');
  }
  return parseResult[1];
};

module.exports = {
  formatTokenForHeader,
  parseToken
};
