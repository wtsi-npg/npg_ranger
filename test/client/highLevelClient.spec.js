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

  it('Success with DNANexus - 1k genomes', ( done ) => {
    var req = new RangerRequest();
    // DNANexus
    var url = 'http://htsnexus.rnd.dnanex.us/v1/reads/1000genomes_low_coverage/NA20276?referenceName=22&start=16100000&end=16105000&format=BAM';

    var expected = {
      bytes: 279561,
      md5: 'd7924c94113ca4e76d64d321f04e6766'
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

  it('Success with DNANexus - platinum', ( done ) => {
    var req = new RangerRequest();
    // DNANexus
    var url = 'http://htsnexus.rnd.dnanex.us/v1/reads/platinum/NA12891?referenceName=chr22&start=16100000&end=16105000&format=BAM';

    var expected = {
      bytes: 814155,
      md5: 'ca298fbb2346bd2df88e647d2e739827'
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
  }, 20000);

  it('Not success with wrong path', ( done ) => {
    var req = new RangerRequest();
    // Google
    var url = 'http://104.196.18.135/readgroupset/CMvnhpKTFhD04eLE-q2yxnU?referenceName=1&start=167856&end=173507&format=BAM';

    req.open('GET', url);

    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.status).toBe(404, 'Correct status code');
        done();
      }
    };

    req.send('');
  }, 5000);

  it('Not success with wrong server', ( done ) => {
    var req = new RangerRequest();
    // Google
    var url = 'http://127.0.0./someData';

    req.open('GET', url);

    req.onreadystatechange = () => {
      if ( req.readyState === 4 ) {
        expect(req.status).toBe(424, 'Correct status code for server address error');
        expect(req.statusMessage).toMatch(/ENOTFOUND/, 'Error message includes statusMessage');
        expect(req.statusMessage).toMatch(/127\.0\.0\./, 'Error message includes the url requested');
        done();
      }
    };

    req.send('');
  }, 5000);
});
