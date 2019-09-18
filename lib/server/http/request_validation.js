"use strict";

/**
* Checks if a duplicate value is passed through in the object.
* Differentiates a duplicate attribute, which would be an array of
* values and should be refused, from multiple regions, which
* would be an array of objects and should be passed through.
* E.g. {"regions":[{ "referenceName" : "chr1" },{ "referenceName" : "chr2"}]}
* would be valid. (Returns undefined.)
* {"format":["bam","cram"]} would be invalid. (Returns 'format'.)
* @return - undefined if valid query, or the failing
* query parameters if invalid.
* @param {object} query - POST query parameters object.
*/

let duplicateAttr = (query) => {
  let isObject = (currentValue) => {
    return typeof currentValue === 'object';
  };
  let dup = Object.keys(query).find(
    (element) => {
      return (isObject(query[element]) &&
              !(Array.isArray(query[element]) && query[element].every(isObject)));
    }
  );
  return dup;
};

module.exports = {
  duplicateAttr
};