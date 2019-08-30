"use strict";

/**
* Checks if a duplicate value is passed through.
* Returns either undefined if the query is valid,
* or the failing query parameters if query is invalid.
* @param query - POST query parameters object
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