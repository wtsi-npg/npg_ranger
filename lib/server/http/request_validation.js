"use strict";

/**
* Checks if a duplicate value is passed through.
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