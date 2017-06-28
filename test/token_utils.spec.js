/* globals describe, it, expect */

"use strict";

let tokenUtils = require('../lib/token_utils');

describe('Errors reported with wrong parameters', () => {
  let bad_params = [
    null,
    undefined,
    1
  ];
  bad_params.forEach( item => {
    it(`throws with wrong parameter type '${item}' `, () => {
      expect(() => {
        tokenUtils.parseToken(item);
      }).toThrowError(/String parameter is required/i);
    });
  });

  bad_params.forEach( item => {
    it(`throws with wrong parameter type '${item}' `, () => {
      expect(() => {
        tokenUtils.formatTokenForHeader(item);
      }).toThrowError(/String parameter is required/i);
    });
  });

  [
    '',
    '    '
  ].forEach( item => {
    it('throws with empty string', () => {
      expect(() => {
        tokenUtils.formatTokenForHeader(item);
      }).toThrowError(/token length must be greater than 0/i);
    });
  });
});

describe('Errors reported when problems with token string', () => {
  it('throws when token string to parse is too long',  () => {
    let long_string = 'Bearer ' + 'a'.repeat(1e3);
    expect( () => {
      tokenUtils.parseToken(long_string);
    }).toThrowError(/is too long/i);
  });
});

describe('Errors are reported when impossible to parse token', () => {
  let bad_formats = [
    '  ',
    ' xxx xxx ',
    'x Bearer xxxx',
    'Bearer',
    ' Bearer ',
    'Bearer ',
    ' Bearer xxxx Bearer xxxx ',
    ' Bearer xxxx xxxx ',
  ];
  bad_formats.forEach(item => {
    it(`throws when malformed token bearer string provided '${item}'`, () => {
      expect(() => {
        tokenUtils.parseToken(item);
      }).toThrowError(/Unexpected format in authorization string/i);
    });
  });
});

describe('Formatted tokens looks like expected', () => {
  let token = 'AAAABBBBCCCC';
  let tokenString = ` Bearer ${token} `;
  it('formats token correctly', () => {
    let formated = tokenUtils.formatTokenForHeader(token);
    expect(formated).toBe(tokenString.trim());
  });
  it('Can parse correctly formated token string', () => {
    let parsed = tokenUtils.parseToken(tokenString);
    expect(parsed).toBe(token);
  });
});
