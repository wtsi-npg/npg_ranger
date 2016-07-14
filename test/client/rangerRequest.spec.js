/* globals describe, expect, it */

"use strict";

const RangerRequest = require('../../lib/client/rangerRequest');

describe('Testing ranger request', () => {
  it('open parameters', ( done ) => {
    var req = new RangerRequest();

    expect(() => {req.open();}).toThrowError(/method is required/);

    var url = 'http://localhost:80/someData';

    expect(() => {req.open('GET');}).toThrowError('url is required');
    expect(() => {req.open('GET', url, 1);}).toThrowError(/async can only be boolean/);
    expect(() => {req.open('GET', url, false);}).toThrowError(/Only async requests are supported/);

    done();
  });

  it('state validation before send', ( done ) => {
    var req = new RangerRequest();
    expect(() => {req.send();}).toThrowError(/The object state must be OPENED/);
    done();
  });
});
