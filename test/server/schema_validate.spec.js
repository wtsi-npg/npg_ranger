/* globals describe, it, expect, beforeAll, afterAll*/

"use strict";

const schemaValid = require('../../lib/server/schema_validate.js');

describe('Valid POST queries', function() {

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
});

describe('Invalid POST parameters - format and accession', function() {

  beforeAll( () => {});
  afterAll( () => {});


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
    it(`testing string definitions for format: ${optionsList.format} and accession: ${optionsList.accession}` , function() {
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
});


describe('POST parameters - regions array', function() {

  beforeAll( () => {});
  afterAll( () => {});
  
  let regionsArray = [ // May remove accession later
    {"regions" : [{"referenceName" : ["ch1"], "start" : 5, "end" : 96}]},
    {"regions" : [{"referenceName" : "ch1", "start" : 96, "end" : 5}]},
    {"regions" : [{"referenceName" : "chr1", "start" : true, "end" : false}]},
    {"regions" : [{"referenceName" : 123, "start" : 5, "end" : 96}]},
    {"regions" : [{"referenceName" : "chr1", "start" : null, "end" : null}]}
  ];

  regionsArray.forEach( optionsList => {
    it(`Invalid queries, testing regions parameters definitions for regions: ${JSON.stringify(optionsList.regions)}` , function() {
      let query = {
        "regions" : optionsList.regions
      };
      expect(() => {schemaValid.validate(query);}).toThrowError();
    });
  });

  it('Query merging - region start and end undefined 1', function() {
    let query = {"regions" : [
      { "referenceName" : "chr2"},
      { "referenceName" : "chr2", "start" : 150, "end" : 200}] };
    let expected = {"regions" : [{ "referenceName" : "chr2"}] };
    expect(schemaValid.validate(query)).toEqual(expected);
  });

  it('Query merging - region start and end undefined 2', function() {
    let query = {"regions" : [
      { "referenceName" : "chr2", "start" : 50, "end" : 100 },
      { "referenceName" : "chr2"}] };
    let expected = {"regions" : [{ "referenceName" : "chr2"}] };
    expect(schemaValid.validate(query)).toEqual(expected);
  });

  it('Query merging - overlapping on both start and end 1', function() {
    let query = {"regions" : [
      { "referenceName" : "chr2", "start" : 50, "end" : 100 },
      { "referenceName" : "chr2", "start" : 25, "end" : 140}] };
    let expected = {"regions" : [{ "referenceName" : "chr2", "start" : 25, "end" : 140}] };
    expect(schemaValid.validate(query)).toEqual(expected);
  });

  it('Query merging - overlapping on both start and end 2', function() {
    let query = {"regions" : [
      { "referenceName" : "chr2", "start" : 25, "end" : 140},
      { "referenceName" : "chr2", "start" : 50, "end" : 100}]};
    let expected = {"regions" : [{ "referenceName" : "chr2", "start" : 25, "end" : 140}] };
    expect(schemaValid.validate(query)).toEqual(expected);
  });

  it('Query merging - overlapping on start', function() {
    let query = {"regions" : [
      { "referenceName" : "chr2", "start" : 1, "end" : 100 },
      { "referenceName" : "chr2", "start" : 25, "end" : 100}] };
    let expected = {"regions" : [{ "referenceName" : "chr2", "start" : 1, "end" : 100}] };
    expect(schemaValid.validate(query)).toEqual(expected);
  });

  it('Query merging - overlapping on start and undefined end', function() {
    let query = {"regions" : [
      { "referenceName" : "chr2", "start" : 1, "end" : 100 },
      { "referenceName" : "chr2", "start" : 25}] };
    let expected = {"regions" : [{ "referenceName" : "chr2", "start" : 1}] };
    expect(schemaValid.validate(query)).toEqual(expected);
  });

  it('Query merging - overlapping on end', function() {
    let query = {"regions" : [
      { "referenceName" : "chr2", "start" : 50, "end" : 100 },
      { "referenceName" : "chr2", "start" : 50, "end": 150}] };
    let expected = {"regions" : [{ "referenceName" : "chr2", "start" : 50, "end" : 150}] };
    expect(schemaValid.validate(query)).toEqual(expected);
  });


  it('Query merging - several regions', function() {
    let query = {"regions" : [
      { "referenceName" : "chr1", "start" : 50, "end" : 100 },
      { "referenceName" : "chr2", "start" : 5, "end" : 10 },
      { "referenceName" : "chr2", "start" : 20, "end" : 100 },
      { "referenceName" : "chr2", "start" : 50, "end": 150}] };
    let expected = {"regions" : [{ "referenceName" : "chr1", "start" : 50, "end" : 100},
                                 { "referenceName" : "chr2", "start" : 5, "end" : 10},
                                 { "referenceName" : "chr2", "start" : 20, "end" : 150}]};
    expect(schemaValid.validate(query)).toEqual(expected);
  });
  
  it('Query merging - several mixed regions', function() {
    let query = {"regions" : [
      { "referenceName" : "chr1", "start" : 50, "end" : 100 },
      { "referenceName" : "chr2", "start" : 50, "end" : 100 },
      { "referenceName" : "chr2", "start" : 5, "end" : 10 },
      { "referenceName" : "chr3", "start" : 50, "end" : 100 },
      { "referenceName" : "chr2", "start" : 120, "end": 150},
      { "referenceName" : "chr3", "start" : 500, "end" : 1000 }]};
    let expected = {"regions" : [{ "referenceName" : "chr1", "start" : 50, "end" : 100},
                                 { "referenceName" : "chr2", "start" : 5, "end" : 10},
                                 { "referenceName" : "chr2", "start" : 50, "end" : 100},
                                 { "referenceName" : "chr2", "start" : 120, "end" : 150},
                                 { "referenceName" : "chr3", "start" : 50, "end" : 100},
                                 {"referenceName" : "chr3", "start" : 500, "end" : 1000}]};
    expect(schemaValid.validate(query)).toEqual(expected);
  });
  
  it('Query merging - two non-overlapping regions', function() {
    let query = {"regions" : [
      { "referenceName" : "chr2", "start" : 50, "end" : 100 },
      { "referenceName" : "chr2", "start" : 150, "end": 200}] };
    let expected = {"regions" : [{ "referenceName" : "chr2", "start" : 50, "end" : 100},
                                 { "referenceName" : "chr2", "start" : 150, "end" : 200}]}; 
    expect(schemaValid.validate(query)).toEqual(expected);
  });
});