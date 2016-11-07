"use strict";

const assert = require('assert');

const SAMTOOLS    = 'samtools';
const BLAT        = 'blat';
const TWO_BIT_EXT = '.2bit';
const INDEX_EXT   = '.fai';
const refPathParser = new RegExp(
  /^(references[\w\.\-\/]+\/)fasta([\w\.\-\/]+)$/,
  'i'
);
const speciesRegExp = new RegExp(
  /^references\/([^\s\/]+)\/([^\s\/]+)\/[^\s]+$/,
  'i'
);

let getReferenceInfo = ( referencePath ) => {
  assert(referencePath, 'referencePath is required');
  let parseResult = speciesRegExp.exec(referencePath);
  let referenceInfo = {
    speciesName:   '',
    referenceName: ''
  };
  if ( parseResult ) {
    if ( parseResult.length > 2 ) {
      let species       = parseResult[1];
      let referenceName = parseResult[2];
      species       = species.replace(/[^a-zA-Z0-9]+/g, ' ');
      referenceName = referenceName.replace(/[^a-zA-Z0-9]+/g, ' ');
      referenceInfo.speciesName = species.trim();
      referenceInfo.referenceName = referenceName.trim();
      referenceInfo.twoBitPath = build2bitPath(referencePath);
    }
  }
  return referenceInfo;
};

let getRange = ( chrData, chr, rangeSize ) => {
  assert(chrData, 'chrData is required');
  assert(chr, 'chr is required');
  assert(rangeSize, 'rangeSize is required');
  let range = null;
  for ( let i = 0 ; i < chrData.length; i++ ) {
    if ( chrData[i][0] == chr ) {
      range = {};
      let chrLength = Number(chrData[i][1]).valueOf();
      let pStart = Math.floor( chrLength / 2) - Math.floor( rangeSize / 2 );
      let pEnd   = Math.floor( chrLength / 2) + Math.floor( rangeSize / 2 );
      range.from = pStart > 1 ? pStart : 1;
      range.to   = pEnd < chrLength ? pEnd : chrLength - 1;
      break;
    }
  }
  return range;
};

let deducePath = ( refPath, toolName, extension ) => {
  assert(refPath, 'parseResult is required');
  assert(toolName, 'toolName is required');
  assert(extension, 'extension is required');
  refPathParser.exec(''); // To remove possible state
  let parseResult = refPathParser.exec(refPath);
  if ( !parseResult ) {
    return null;
  }
  assert(
    parseResult.length > 2,
    'Result of parsing reference path is shorter than expected'
  );
  assert(
    parseResult[1] && parseResult[2],
    'Result of parsing reference does not match expected components'
  );
  let pre = parseResult[1];
  let suf = parseResult[2];
  let twoBitPath = pre + toolName + suf + extension;
  return twoBitPath;
};

let build2bitPath = ( refPath ) => {
  return deducePath(refPath, BLAT, TWO_BIT_EXT);
};

let buildIndexPath = ( refPath ) => {
  return deducePath(refPath, SAMTOOLS, INDEX_EXT);
};

let parseIndexData = ( indexData ) => {
  assert(indexData, 'indexData is required');
  assert(typeof indexData === 'string', 'indexData must be a string');
  let asLines = indexData.split(/\r?\n/);
  let linesAsColumns = [];
  for ( let i = 0; i < asLines.length; i++ ) {
    if ( asLines[i].trim() !== '') {
      linesAsColumns.push(asLines[i].split(/[\s]+/).slice(0,2));
    }
  }
  return linesAsColumns;
};

module.exports = {
  buildIndexPath: buildIndexPath,
  getRange:       getRange,
  getReferenceInfo: getReferenceInfo,
  parseIndexData: parseIndexData
};
