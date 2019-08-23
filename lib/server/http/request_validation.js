"use strict";

let duplicateAttr = (query) => {
  let isObject = (currentValue) => {
    return typeof currentValue === 'object';
  };
  let testing = Object.keys(query).find(
    (element) => {
      /* console.log("current element is: ");
      console.log(element);
      console.log(query[element]);
      console.log("Is element an object: ");
      console.log(typeof query[element] === 'object'); */
      return (typeof query[element] === 'object' &&
              !(Array.isArray(query[element]) && query[element].every(isObject)));
    } // TODO
  );
  /* console.log("stuff = ");
  console.log(testing);
  console.log(typeof testing === 'undefined'); */
  return testing;
};

module.exports = {
  duplicateAttr
};