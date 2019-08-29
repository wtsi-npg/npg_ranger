"use strict";

const jsonLintSchema  = require('@prantlf/jsonlint/lib/validator');
const htsSchema       = require('./htsget_post.schema.json');
// const ServerHttpError = require('./http/error.js'); // 
// const HttpError       = ServerHttpError.HttpError; // 

let startEndValidator = (start, end) => {
  if (start >= end) {
    
    throw new RangeError('start is greater or equal to end.'); // TODO  Throw InvalidRange instead?
  }
};

let regionMerger = (regions) => {
  let rgnHasher = (regions) => {
    let newHash = {}; // In the form of {"a" : [[1,2],[3,4]]}
    regions.forEach( indvRegion => {
      let refName = indvRegion.referenceName;
      if (typeof newHash[refName] === 'undefined') {
        newHash[refName] = [];
      }
      newHash[refName].push([indvRegion.start, indvRegion.end]);
    });
    console.log('new hash is:');
    console.log(newHash);
    console.log('Sorting the hash:');
    for (let val of Object.values(newHash)){
      val.sort((a,b) => a[0] - b[0]);
    }
    console.log(newHash);
    return newHash;
  };
  let rgnCollapse = (inputRgn) => {
    var outputRgn = {};
    console.log('Entering rngCollapse, with input region of: ');
    console.log(inputRgn);
    for (const [key, value] of Object.entries(inputRgn)) { 
      // Will be getting (chr1 : [[5,10],[7,15]]) styled
      // inputs from the entirety of the input region
      console.log('key and value of the inputRgn');
      console.log(key);
      console.log(value);
      value.forEach ( inputPair => {
        // for a given key, loop through each value
        // so input would be [5,10], then [7,15]
        // all for chr1
        if (typeof outputRgn[key] === 'undefined') {
          outputRgn[key] = [];
          outputRgn[key].push(inputPair);
        } else {
          if (typeof inputPair[0] === 'undefined' || inputPair[0] < outputRgn[key][outputRgn[key].length - 1][1]) { // TODO check if this is clean enough?
            if (typeof inputPair[0] === 'undefined' || 
                (inputPair[0] <  outputRgn[key][outputRgn[key].length - 1][0])) {
              outputRgn[key][outputRgn[key].length - 1][0] = inputPair[0];
            }
            if (typeof inputPair[1] === 'undefined' ||  inputPair[1] > outputRgn[key][outputRgn[key].length - 1][1]) {
              outputRgn[key][outputRgn[key].length - 1][1] = inputPair[1];
            }
            
          } else {
            outputRgn[key].push(inputPair);
          }
        }
      });
    }
    console.log('outputRgn');
    console.log(outputRgn);
    return outputRgn;
  };
  let rgnRebuild = (hash) => {
    var parsedRegion = [];
    for (const [key, value] of Object.entries(hash)) { 
      // let finalObj = {};
      // let refNameObj = {};
      // refName.referenceName = key;
      value.forEach( pair => {
        let tempObj = {};
        tempObj.referenceName = key; // TODO worth the effort to avoid reassigning?
        tempObj.start = pair[0];
        tempObj.end = pair[1];
        parsedRegion.push(tempObj);
      });
    }
    console.log('final response: ');
    console.log(parsedRegion);
    return parsedRegion;
  };
  
  let hashedRegion = rgnCollapse(rgnHasher(regions));
  let outputRegion = rgnRebuild(hashedRegion);
  return outputRegion;
};

let validate = (test) => {
  let POSTSchema = jsonLintSchema.compile(JSON.stringify(htsSchema));
  let tmpTest = new POSTSchema(test);
  let validRegion = [];
  console.log(tmpTest);
  console.log(tmpTest.regions);
  try {
    if (tmpTest.regions) {
      let regionsClone = JSON.parse(JSON.stringify(tmpTest.regions));
      //regionsClone.sort(function (a, b) {
      //  return (a.referenceName).localeCompare(b.referenceName);
      //});

      regionsClone.forEach( regionList => {
        console.log(regionList);
        startEndValidator(regionList.start, regionList.end); // TODO - doublecheck
      });
      validRegion = regionMerger(regionsClone);
    }
  } catch (e) { // TODO fix up this region
    // let err = new HttpError(tmpTest, 400, 'InvalidRange: start cannot be greater or equal to end.');
    // console.log(err);
    // throw err;
    throw e;
  }
  tmpTest.regions = validRegion;
  console.log("validation ends with: ");
  console.log(validRegion);
  console.log(tmpTest);
  return tmpTest;
};

module.exports = {
  validate
};