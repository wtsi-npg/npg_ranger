/* globals describe, expect, it */

"use strict";

const md5 = require('js-md5');

const RangerRequest = require('../../lib/client/rangerRequest');

describe('Testing high level', () => {
  it('Success with Google', ( done ) => {
    var req = new RangerRequest();
    // Google
    var url = 'http://104.196.18.135/readgroupsets/CMvnhpKTFhD04eLE-q2yxnU?referenceName=1&start=167856&end=173507&format=BAM';

    var expected = {
      bytes: 63778,
      md5: '94a3fc74146898aab618cad666a3a54a'
    };

    req.open('GET', url);

    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.status === 200 || req.status === 206).toBe(true);
        expect(req.response.byteLength).toBe(expected.bytes, 'Got correct number of bytes');
        expect(md5(new Uint8Array(req.response))).toBe(expected.md5, 'Data hash matches expected');
        done();
      }
    };

    req.send('');
  }, 5000);

  it('Success with Google, but with different data', ( done ) => {
    var req = new RangerRequest();
    // Google
    var url = 'http://104.196.18.135/readgroupsets/CMvnhpKTFhD04eLE-q2yxnU?referenceName=1&start=160000&end=165000&format=BAM';

    var expected = {
      bytes: 57977,
      md5: '41a7756429c37a78b65b5aa4a5891e54'
    };

    req.open('GET', url);

    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.status === 200 || req.status === 206).toBe(true);
        expect(req.response.byteLength).toBe(expected.bytes, 'Got correct number of bytes');
        expect(md5(new Uint8Array(req.response))).toBe(expected.md5, 'Data hash matches expected');
        done();
      }
    };

    req.send('');
  }, 5000);

  it('Not sucess with wrong path', ( done ) => {
    var req = new RangerRequest();
    // Google
    var url = 'http://104.196.18.135/readgroupset/CMvnhpKTFhD04eLE-q2yxnU?referenceName=1&start=167856&end=173507&format=BAM';

    req.open('GET', url);

    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.status === 200 || req.status === 206).not.toBe(true);
        done();
      }
    };

    req.send('');
  }, 5000);
});