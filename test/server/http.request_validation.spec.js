/* globals describe, it, expect, beforeAll, afterAll*/

"use strict";

const reqValid = require('../../lib/server/http/request_validation.js');

describe('Temp - all tests', function() {

  beforeAll( () => {});
  afterAll( () => {});

  it('Pass valid request, no objects', function() {
    let tmp = reqValid.duplicateAttr({"format" : "bam"});
    expect(typeof tmp !== 'undefined').toBe(false);
  });

  it('Invalid array, does not contain objects', function() {
    let tmp = reqValid.duplicateAttr(({
      "format" : "bam",
      "regions" : ["foo", "bar"]
    }));
    expect(typeof tmp !== 'undefined').toBe(true);
  });

  it('Pass valid several queries, no objects', function() {
    let tmp = reqValid.duplicateAttr({"format" : "bam", "tags" : "RG"});
    expect(typeof tmp !== 'undefined').toBe(false);
  });

  it('Valid array, holds objects', function() {
    let tmp = reqValid.duplicateAttr({
      "format" : "bam",
      "regions" : [{"referenceName" : "chr1"}]
    });
    expect(typeof tmp !== 'undefined').toBe(false);
  });

  it('Object passed in request, not array', function() {
    let tmp = reqValid.duplicateAttr({
      "format" : "bam",
      "regions" : {"referenceName" : "chr1"}
    });
    expect(typeof tmp !== 'undefined').toBe(true);
  });
});

