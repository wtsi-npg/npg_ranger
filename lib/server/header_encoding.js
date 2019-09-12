"use strict";

const zlib = require('zlib');
const fs = require('fs');
const buffer = require('buffer');

//const gzip = zlib.createGzip();
// const inp = fs.createReadStream('fileHeaderInput.txt'); // 15962
//const outGz = fs.createWriteStream('fileHeaderInput.txt.gz'); // 2114
//const outB64 = fs.createWriteStream('fileHeaderInputB64.txt'); // 
//inp.pipe(gzip).pipe(out);

let genTestData = async () => {
  return new Promise((resolve, reject) => {
    // let fullRegions = {"regions" : []};
    let tempRegArray = [];
    let makeRegions = ( regArray ) => {
      let temp;
      for (let i = 0; i < 1000; i++) {
        temp = {};
        temp.referenceName = "chr" + i;
        if (i % 4 == 0) {
          temp.start = Math.floor(Math.random() * 10000) + 1;
        }
        if (i % 6 == 0) {
          temp.end = Math.floor(Math.random() * 10000) + 1;
        }
        regArray.push(temp);
      }
    };
    makeRegions( tempRegArray );
    // console.log( tempRegArray );
    let output = JSON.stringify(tempRegArray);
    console.log(output);
    // fs.writeFile("fileHeaderInput.txt", output, function(err) {
    //   if (err) {
    //     console.log('aa');
    //   }
    // });
    resolve(new Buffer.from(output, 'utf8'));
  });
};


let readFile = async (filename) => { //into buffer
  return new Promise((resolve, reject) => {
    fs.readFile(filename, function(err, data) {
      if (err) { 
        console.log('aaa');
        reject(err);
      } else {
        //console.log(JSON.parse(data));
        resolve(new Buffer.from(data, 'utf8'));
      }
    });
  });
};

let gzipFile = async (inputBuffer) => {
  return new Promise((resolve, reject) => {
    zlib.gzip(inputBuffer, function (error, result) {
      if (error) { throw error; }
      let value = result;
      //fs.writeFile("fileHeaderInputGz.txt", value, (err) => {
      //  if (err) { throw err; }
      //});
      //console.log(inputBuffer.toString());
      //console.log(temp);
      //console.log(value);
      resolve(value);
    });
  });
};

let unzipFile = async (inputBuffer) => { 
  return new Promise((resolve, reject) => {
    zlib.unzip(inputBuffer, function (error, result) {
      if (error) { throw error; }
      let value = result;
      fs.writeFile("fileHeaderInputGzToTxt.txt", value, (err) => {
        if (err) { throw err; }
      });
      //console.log(value);
      resolve(value);
    });
  });
};


let b64ConvFile = async (inputBuffer) => {
  return new Promise((resolve, reject) => {
    let base64Input = new Buffer.from(inputBuffer, 'binary').toString('base64');
    //let base64Input = buffer.transcode(inputBuffer, 'binary', 'base64');
    fs.writeFile("fileHeaderInputB64.txt", base64Input, (err) => {
      if (err) { throw err; }
    });
    console.log('thing');
    console.log();
    resolve(base64Input);
  });
};

let b64DeconvFile = async (inputBuffer) => {
  return new Promise((resolve, reject) => {
    let asciiInput = new Buffer.from(inputBuffer, 'base64');
    fs.writeFile("fileHeaderInputB64toAscii.txt", asciiInput, (err) => {
      if (err) { throw err; }
    });
    console.log('thing2');
    console.log(asciiInput);
    resolve(asciiInput);
  });
};

let executeStuff = async () => {
  let buff = await readFile('fileHeaderInput.txt'); //Read input file, turn into a buffer.
  let zip = await gzipFile(buff); //Gzip the buffer.
  let unzip = await unzipFile(zip); //unzip zipped buffers
  let b64 = await b64ConvFile(zip); //Convert zipped buffer into b64
  let ascii = await b64DeconvFile(b64); //Deconvert b64 zipped buffer into normal buffer
  let temp = await unzipFile(ascii);
  //Something stupid is going on here - due to the gzip/outGZ stuff, the file is wiped.;
  //await b64ConvFile(); //So the file is wiped first, then b64 tries to convert it, fails, and then gzip remakes it.
  return Promise.resolve();
};

executeStuff();

module.exports = {
  
};