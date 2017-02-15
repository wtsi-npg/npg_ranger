/* globals describe, it, expect */

"use strict";

const uriUtils = require('../../lib/client/uriUtils');

describe('Testing processing data URIs', () => {
  let expected = 'some data';
  let encoded = new Buffer(expected, 'ascii').toString('base64');
  let dataURI = 'data:text/plain;charset=utf-8;base64,' + encoded;

  it('can get buffer from base64 data URI', () => {
    let buffer = uriUtils.procDataURI(dataURI);
    expect(typeof buffer).toEqual('object');
    expect(buffer instanceof Buffer).toBe(true);
    expect(buffer.toString()).toEqual(expected);
  });
});

describe('Fails with missing or unknown encoding', () => {
  let wrongURIs = [
    'data:text/plain;charset=utf-8;unknown,xxxx',
    'data:text/plain;charset=utf-8,xxxx'
  ];

  wrongURIs.forEach( uri => {
    it(`fails with "${uri}"`, () => {
      expect(()=> {
        uriUtils.procDataURI(uri);
      }).toThrowError(/(Unknown encoding)|("encoding" must be a valid string)/);
    });
  });
});

describe('Fails with non-data URIs', () => {
  let wrongURIs = [
    'file://myfile.txt',
    'http://localhost/myfile.txt',
    'https://localhost/myfile.txt',
    'random text'
  ];

  wrongURIs.forEach( uri => {
    it(`fails with "${uri}"`, () => {
      expect(()=> {
        uriUtils.procDataURI(uri);
      }).toThrowError(/Unable to decode uri/);
    });
  });
});
