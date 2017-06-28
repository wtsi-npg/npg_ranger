"use strict";

const assert = require('assert');
const util   = require('util');

const TOKEN_REGEX          = /^\s*Bearer\s+(\S+)\s*$/i;
const TOKEN_BEARER_FORMAT  = `Bearer %s`;
const TOKEN_BEARER_MAX_LEN = 1000;

let formatTokenForHeader = token => {
  assert(typeof token === 'string', 'String parameter is required');
  token = token.trim();
  assert(token.length > 0, 'token length must be greater than 0');
  return util.format(TOKEN_BEARER_FORMAT, token);
};

let parseToken = string => {
  assert(typeof string === 'string', 'String parameter is required');
  if ( string.length > 1000 ) {
    throw new Error(
      `String provided for parsing token is too long (>${TOKEN_BEARER_MAX_LEN})`
    );
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
