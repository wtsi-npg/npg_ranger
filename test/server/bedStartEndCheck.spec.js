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

  it('No end given, TEMP: ASSUME LARGE NUMBER, FIND VALID END REGION IN FUTURE', function() { // TODO modify test to be correct for end values
    let region = [{"referenceName": "test", "start": 100}];
    startEndConstructor(region);
    expect(region).toEqual([{"referenceName": "test", "start": 100, "end": 30000}]);
  });

  it('Neither start or end given, assume for both', function() { // TODO modify test to be correct for end values
    let region = [{"referenceName": "test"}];
    startEndConstructor(region);
    expect(region).toEqual([{"referenceName": "test", "start": 0, "end": 30000}]);
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