/*globals describe, it, expect, beforeAll, afterAll*/

"use strict";

const headerEncode = require('../../lib/server/header_encoding.js');

describe('Temp Describe Name', function() {

  beforeAll( () => {
    //let testData = headerEncode.genTestData();
  });
  afterAll( () => {});

  let genTestData = async (refCount, rangeSize) => { // generate test data and return it as a buffer.
    return new Promise(( resolve ) => {
      // let fullRegions = {"regions" : []};
      let tempRegArray = [];
      let makeRegions = ( regArray ) => {
        let temp;
        for (let i = 0; i < refCount; i++) {
          temp = {};
          temp.referenceName = "chr" + i;
          if (i % 4 == 0) {
            temp.start = Math.floor(Math.random() * rangeSize) + 1;
          }
          if (i % 6 == 0) {
            temp.end = Math.floor(Math.random() * rangeSize) + 1;
          }
          regArray.push( temp );
        }
      };
      makeRegions( tempRegArray );
      let output = JSON.stringify( tempRegArray );
      resolve(new Buffer.from(output, 'utf8'));
    });
  };

  it('Running on randomly generated large data', async function() {
    let testData = await genTestData(2000, 1000000);
    let convertedData = await headerEncode.fullEncoding(testData);
    let deconvertedData = await headerEncode.fullDecoding(convertedData);
    expect(testData.equals(deconvertedData)).toBe(true);
  });

  it('Running on randomly generated small data', async function() {
    let testData = await genTestData(10, 1000000);
    let convertedData = await headerEncode.fullEncoding(testData);
    let deconvertedData = await headerEncode.fullDecoding(convertedData);
    expect(testData.equals(deconvertedData)).toBe(true);
  });

  it('Empty array input', async function() {
    let testData = new Buffer.from([]);
    let convertedData = await headerEncode.fullEncoding(testData);
    let deconvertedData = await headerEncode.fullDecoding(convertedData);
    expect(testData.equals(deconvertedData)).toBe(true);
  });
  
  it('Empty array input', async function() {
    let testData = new Buffer.from([]);
    let str = 'H4sIAAAAAAAAA4uuVipKTUstSs1LTvVLzE1VslJKzigyVKrVwSphpKSjVFySWFSiZGVpaamjlJqXomRlaGBgQFi9EVAVVIMRUEdtLAAiqs+mewAAAA=='
    let temp = new Buffer.from(str);
    let deconvertedData = await headerEncode.fullDecoding(str);
    console.log(deconvertedData);
    expect(testData.equals(deconvertedData)).toBe(true);
  });
 
});