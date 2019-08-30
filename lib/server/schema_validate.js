"use strict";

const jsonLintSchema  = require('@prantlf/jsonlint/lib/validator');
const htsSchema       = require('./htsget_post.schema.json');

const POSTSchema = jsonLintSchema.compile(JSON.stringify(htsSchema));

let _startEndValidator = (start, end) => {
  if (start >= end) {
    throw new Error('Range end should be bigger than start');
  }
};

let _regionMerger = (regions) => {
  let rgnHasher = (regions) => {
    let newHash = {}; // In the form of {"a" : [[1,2],[3,4]]}
    regions.forEach( indvRegion => {
      let refName = indvRegion.referenceName;
      if (typeof newHash[refName] === 'undefined') {
        newHash[refName] = [];
      }
      newHash[refName].push([indvRegion.start, indvRegion.end]);
    });
    for (let val of Object.values(newHash)) {
      val.sort((a,b) => a[0] - b[0]);
    }
    return newHash;
  };
  let rgnCollapse = (inputRgn) => {
    let outputRgn = {};
    /* jshint -W083 */
    for (let [key, value] of Object.entries(inputRgn)) {
      value.forEach ( inputPair => {
        if (typeof outputRgn[key] === 'undefined') {
          outputRgn[key] = [];
          outputRgn[key].push(inputPair);
        } else {
          // Merges the regions if fits the criteria
          if (typeof outputRgn[key][outputRgn[key].length - 1][1] === 'undefined' ||
              typeof inputPair[0] === 'undefined' ||
              inputPair[0] < outputRgn[key][outputRgn[key].length - 1][1]) {
            if (typeof inputPair[0] === 'undefined' ||
                inputPair[0] <  outputRgn[key][outputRgn[key].length - 1][0]) {
              outputRgn[key][outputRgn[key].length - 1][0] = inputPair[0];
            }
            if (typeof inputPair[1] === 'undefined' ||
                inputPair[1] > outputRgn[key][outputRgn[key].length - 1][1]) {
              outputRgn[key][outputRgn[key].length - 1][1] = inputPair[1];
            }
          } else {
            outputRgn[key].push(inputPair);
          }
        }
      });
    }
    /* jshint +W083 */
    return outputRgn;
  };
  let rgnRebuild = (hash) => {
    let parsedRegion = [];
    /* jshint -W083 */
    for (let [key, value] of Object.entries(hash)) {
      value.forEach( pair => {
        let tempObj = {};
        tempObj.referenceName = key;
        // Do not include 'undefined' for the output
        if (pair[0]) {tempObj.start = pair[0];}
        if (pair[1]) {tempObj.end = pair[1];}
        parsedRegion.push(tempObj);
      });
    }
    /* jshint +W083 */
    return parsedRegion;
  };
  let hashedRegion = rgnCollapse(rgnHasher(regions));
  let outputRegion = rgnRebuild(hashedRegion);
  return outputRegion;
};

/**
* Validate that the POST query passed through is Valid.
* First checks if input fits the htsget schema file,
* then checks and merges, if required,
* the region parameter if present.
* Return the valid (and if suitable, merged) POST query.
* @param POSTQuery - a valid JSON format POST query.
*/
let validate = (POSTQuery) => {
  let validQuery;
  let validRegion;
  try {
    validQuery = new POSTSchema(POSTQuery);
    validRegion = [];
  } catch (e) {
    throw e;
  }
  try {
    if (validQuery.regions) {
      let regionsClone = JSON.parse(JSON.stringify(validQuery.regions));
      regionsClone.sort = (a, b) => {
        return (a.referenceName).localeCompare(b.referenceName);
      };
      regionsClone.forEach( regionList => {
        _startEndValidator(regionList.start, regionList.end);
      });
      validRegion = _regionMerger(regionsClone);
      validQuery.regions = validRegion;
    }
  } catch (e) {
    throw e;
  }
  return validQuery;
};

module.exports = {
  validate
};