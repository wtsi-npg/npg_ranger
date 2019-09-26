"use strict";

const zlib = require('zlib');

let _gzipFile = async (inputBuffer) => {
  return new Promise((resolve, reject) => {
    zlib.gzip(inputBuffer, function( err, result ) {
      if ( err ) { reject( err ); }
      let value = result;
      resolve( value );
    });
  });
};

let _unzipFile = async (inputBuffer) => {
  return new Promise((resolve, reject) => {
    zlib.unzip(inputBuffer, function( err, result ) {
      if ( err ) { reject( err ); }
      let value = result;
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

module.exports = {
  fullEncoding: fullEncoding,
  fullDecoding: fullDecoding
};