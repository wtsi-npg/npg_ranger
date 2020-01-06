/* globals describe, it, expect */

"use strict";

const bedStrConstructor = require('../../lib/server/model.js').buildBedFile;


describe('Test the string constructor for .bed files', function() {
  
  it('Valid region single element array', function() {
    let region = [{"referenceName": "test", "start": 200, "end": 600}];
    let builtStr = bedStrConstructor(region);
    expect(builtStr).toBe("test\t200\t600\n");
  });
  
  it('Valid region several element array', function() {
    let region = [
      {"referenceName": "test", "start": 200, "end": 600},
      {"referenceName": "abc", "start": 1500, "end": 3200},
      {"referenceName": "phix", "start": 50, "end": 400}
                 ];
    let builtStr = bedStrConstructor(region);
    expect(builtStr).toBe("test\t200\t600\nabc\t1500\t3200\nphix\t50\t400\n");
  });

  // Maybe the StrConstructor should put the inputted array/string through a linter/ validator?

  it('Invalid empty input', function() {
    let region;
    expect( () => { bedStrConstructor(region); }).toThrowError();
  });

  it('Invalid non-array input', function() {
    let region = "referenceName:test";
    expect( () => { bedStrConstructor(region); }).toThrowError();
  });

});