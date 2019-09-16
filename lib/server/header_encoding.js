"use strict";

const zlib = require('zlib');

/*
const fs = require('fs');

let readFile = async (filename) => { // into buffer
  return new Promise(( resolve, reject ) => {
    fs.readFile(filename, function(err, data) {
      if ( err ) {
        // console.log('error on readfile');
        reject( err );
      } else {
        //console.log(JSON.parse(data));
        resolve(new Buffer.from(data, 'utf8'));
      }
    });
  });
};
*/

let _gzipFile = async (inputBuffer) => {
  return new Promise((resolve, reject) => {
    zlib.gzip(inputBuffer, function( err, result ) {
      if ( err ) { reject( err ); }
      let value = result;
      // fs.writeFile("fileHeaderInputGz.txt", value, (err) => {
      //   if (err) { throw err; }
      // });
      // console.log(inputBuffer.toString());
      // console.log(temp);
      // console.log(value);
      resolve( value );
    });
  });
};

let _unzipFile = async (inputBuffer) => {
  return new Promise((resolve, reject) => {
    zlib.unzip(inputBuffer, function( err, result ) {
      if ( err ) { reject( err ); }
      let value = result;
      // fs.writeFile("fileHeaderInputGzToTxt.txt", value, (err) => {
      //   if (err) { throw err; }
      // });
      // console.log(value);
      resolve( value );
    });
  });
};


let _b64ConvFile = async (inputBuffer) => {
  return new Promise((resolve, reject) => {
    try {
      let base64Input = new Buffer.from(inputBuffer, 'binary').toString('base64');
      resolve( base64Input );
    } catch ( err ) {
      reject( err );
    }
    // let base64Input = buffer.transcode(inputBuffer, 'binary', 'base64');
    // fs.writeFile("fileHeaderInputB64.txt", base64Input, (err) => {
    //   if (err) { throw err; }
    // });
  });
};

let _b64DeconvFile = async (inputBuffer) => {
  return new Promise((resolve, reject) => {
    try {
      let asciiInput = new Buffer.from(inputBuffer, 'base64');
      resolve( asciiInput );
    } catch ( err ) {
      reject( err );
    }
    // fs.writeFile("fileHeaderInputB64toAscii.txt", asciiInput, (err) => {
    //   if (err) { throw err; }
    // });
    // console.log(asciiInput);
  });
};

let fullEncoding = async (bufferedData) => {
  let zippedData = await _gzipFile(bufferedData); // Gzip the buffer.
  let b64EncData = await _b64ConvFile(zippedData); // Convert zipped buffer into b64
  return b64EncData; // Return zipped, b64 encoded buffer.
};

let fullDecoding = async (encodedBuffer) => {
  let zippedData = await _b64DeconvFile(encodedBuffer); // Deconvert b64 zipped buffer into normal zipped buffer
  let originalData = await _unzipFile(zippedData);
  return originalData; // Returns the original string in a buffer.
};

/*
let temp = async () => {
  let buff = await genTestData();
  console.log(buff);
  let val = await fullEncoding(buff);
  let val2 = await fullDecoding(val);
  return buff.equals(val2);
};
temp();
*/

module.exports = {
  fullEncoding: fullEncoding,
  fullDecoding: fullDecoding
};