"use strict";

const jsonLintSchema  = require('@prantlf/jsonlint/lib/validator');
const htsSchema       = require('./htsget_post.schema.json');

let startEndValidator = (start, end) => {
  if (start >= end) {
    throw new InvalidRange('start is greater or equal to end.');
  }
};

let validate = (test) => {
  let POSTSchema = jsonLintSchema.compile(JSON.stringify(htsSchema));
  let tmpTest = new POSTSchema(test);
  console.log(tmpTest.regions);
  try {
    if (tmpTest.regions){
      tmpTest.regions.forEach( regionList => {
        startEndValidator(regionList.start, regionList.end)
      });
    }
  } catch (e) { // fix up this region
    console.log(e.toString());
    throw e;
  }

  return tmpTest;
};

module.exports = {
  validate
};