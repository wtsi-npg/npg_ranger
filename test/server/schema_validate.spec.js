/* globals describe, it, expect, beforeAll, afterAll*/

"use strict";

const schemaValid = require('../../lib/server/schema_validate.js');

describe('Temp', function() {

  beforeAll( () => {});
  afterAll( () => {});

  it('confirm spec example', function() {
    let query = {
      "format" : "bam",
      "fields" : ["QNAME", "FLAG", "RNAME", "POS", "CIGAR", "SEQ"],
      "tags" : ["RG"],
      "notags" : ["OQ"],
      "regions" : [
        { "referenceName" : "chr1" }
      ]
    };
    expect(schemaValid.validate(query)).toBe(query);
  });


  it('similar valid query', function() {
    let query = {
      "format" : "cram",
      "fields" : ["RNAME", "POS", "CIGAR", "SEQ"],
      "tags" : ["RG"],
      "notags" : ["OQ"],
      "regions" : [
        { "referenceName" : "chr3", "start" : 999, "end" : 1023 }
      ]
    };
    expect(schemaValid.validate(query)).toBe(query);
  });

  it('valid when query is empty', function() {
    let query = {};
    expect(schemaValid.validate(query)).toBe(query);
  });

  let stringArray = [ // May remove accession later
    {format : 123, accession : "XYZ123456"},
    {format : "bam", accession : 123456789},
    {format : ["bam"], accession : "XYZ123456"},
    {format : "bam", accession : ["XYZ123456"]},
    {format : null, accession : "XYZ123456"},
    {format : "bam", accession : false}
    // {format : {"bam"}, accession : {"XYZ123456"}}
  ];

  stringArray.forEach( optionsList => {
    it(`testing string definitions for format:${optionsList.format} and accession:${optionsList.accession}` , function() {
      let query = {
        "format" : optionsList.format,
        "accession" : optionsList.accession
      };
      expect(() => {schemaValid.validate(query);}).toThrowError();
    });
  });

  let arraysArray = [
    {"fields" : "QNAME",
     "tags" : ["RG"],
     "notags" : ["OQ"]},
    {"fields" : ["QNAME", "FLAG", "RNAME", "POS", "CIGAR", "SEQ"],
     "tags" : 123,
     "notags" : ["OQ"]},
    {"fields" : ["QNAME", "FLAG", "RNAME", "POS", "CIGAR", "SEQ"],
     "tags" : ["RG"],
     "notags" : true},
    {"fields" : {"QNAME" : 1},
     "tags" : {"RG" : 2},
     "notags" : {"QQ" : 3}},
    {"fields" : ["FAKE", "WRONG"],
     "tags" : ["RG"],
     "notags" : ["OQ"]},
  ];

  arraysArray.forEach( optionsList => {
    it(`testing array definitions for fields: ${JSON.stringify(optionsList.fields)}, 
       tags: ${JSON.stringify(optionsList.tags)} 
       and notags: ${JSON.stringify(optionsList.tags)}`, function() {
      let query = {
        "fields" : optionsList.fields,
        "tags" : optionsList.tags,
        "notags" : optionsList.notags
      };
      expect(() => {schemaValid.validate(query);}).toThrowError();
    });
  });

  
  let regionsArray = [ // May remove accession later
    {"regions" : [{"referenceName" : ["ch1"], "start" : 5, "end" : 96}]},
    {"regions" : [{"referenceName" : "ch1", "start" : 96, "end" : 5}]},
    {"regions" : [{"referenceName" : "chr1", "start" : true, "end" : false}]},
    {"regions" : [{"referenceName" : 123, "start" : 5, "end" : 96}]},
    {"regions" : [{"referenceName" : "chr1", "start" : null, "end" : null}]}
  ];

  regionsArray.forEach( optionsList => {
    it(`testing regions parameters definitions for regions:${JSON.stringify(optionsList.regions)}` , function() {
      let query = {
        "regions" : optionsList.regions
      };
      expect(() => {schemaValid.validate(query);}).toThrowError();
    });
  });


});