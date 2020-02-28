/* globals describe, it, expect */

"use strict";

const startEndConstructor = require('../../lib/server/model.js').checkStartEnd;


describe('Test the .bed file start and end regions', function() {
  
  it('No start given, turn to 0', function() {
    let region = [{"referenceName": "test", "end": 600}];
    startEndConstructor(region);
    expect(region).toEqual([{"referenceName": "test", "start": 0, "end": 600}]);
  });  

  it('Valid start and end given, no effect', function() {
    let region = [{"referenceName": "test", "start": 100, "end": 500}];
    startEndConstructor(region);
    expect(region).toEqual([{"referenceName": "test", "start": 100, "end": 500}]);
  });

  it('Valid start and end given 2, no effect', function() {
    let region = [{"referenceName": "test", "start": 1250, "end": 5000}];
    startEndConstructor(region);
    expect(region).toEqual([{"referenceName": "test", "start": 1250, "end": 5000}]);
  });

  it('No end given, do not modify', function() {
    let region = [{"referenceName": "test", "start": 100}];
    startEndConstructor(region);
    expect(region).toEqual([{"referenceName": "test", "start": 100}]);
    expect("end" in region).toBe(false);
  });

  it('Neither start or end given, assume for start, do not modify end', function() {
    let region = [{"referenceName": "test"}];
    startEndConstructor(region);
    expect(region).toEqual([{"referenceName": "test", "start": 0}]);
    expect("end" in region).toBe(false);
  });

  it('Invalid empty input', function() {
    let region;
    expect( () => { startEndConstructor(region); }).toThrowError();
  });

  it('Invalid non-array input', function() {
    let region = "referenceName:test";
    expect( () => { startEndConstructor(region); }).toThrowError();
  });

});