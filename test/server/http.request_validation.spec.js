/* globals describe, it, expect, beforeAll, afterAll*/

"use strict";

const reqValid = require('../../lib/server/http/request_validation.js');

describe('POST request parameters testing', function() {

  beforeAll( () => {});
  afterAll( () => {});

  it('Pass valid request, no objects', function() {
    let tmp = reqValid.duplicateAttr({"format" : "bam"});
    expect(typeof tmp !== 'undefined').toBe(false);
  });

  it('Pass empty request', function() {
    let tmp = reqValid.duplicateAttr({});
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

  it('Valid array, holds basic object', function() {
    let tmp = reqValid.duplicateAttr({
      "format" : "bam",
      "regions" : [{"referenceName" : "chr1"}]
    });
    expect(typeof tmp !== 'undefined').toBe(false);
  });

  it('Valid array, holds more complex object', function() {
    let tmp = reqValid.duplicateAttr({
      "format" : "bam",
      "regions" : [{"referenceName" : "chr2", "start" : 3, "end" : 10}]
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

