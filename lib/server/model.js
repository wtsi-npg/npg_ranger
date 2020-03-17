"use strict";

/**
 * model module.
 * @module server/model
 *
 * @description Model for the streaming HTTP server.
 *
 * @example <caption>Example usage of the model.</caption>
 *   const RangerModel = require('../lib/server/model.js');
 *   // Create a new object.
 *   let m = new RangerModel();
 *   // process data
 *   m.process(query, destination, endResponseClback);
 *
 * @copyright Genome Research Limited 2017
 */

const assert   = require('assert');
const child    = require('child_process');
const fse      = require('fs-extra');
const path     = require('path');
const LOGGER   = require('winston');
const pipeline = require('./pipeline.js');
const config   = require('../config.js');

const SAMTOOLS_COMMAND     = 'samtools';
const FREEBAYES_COMMAND    = 'freebayes';
const DEFAULT_FORMAT       = 'BAM';

const TEMPDIR_FILEPATH     = '/bedFiles';
const BED_FILETYPE       = '.bed';

let _cleanup = ( dir ) => {
  fse.remove(dir, ( err ) => {
    if ( err ) {
      LOGGER.warn(`Failed to remove ${dir}: ${err}`);
    }
  });
};

let buildBedFile = (regions) => {
  assert(regions);
  assert(Array.isArray(regions));
  let bedString = '';
  regions.forEach ( region => {
    bedString += region.referenceName + "\t" + region.start + "\t" + region.end + "\n";
  });
  return bedString;
};

let checkStartEnd = (regions) => {
  assert(regions);
  assert(Array.isArray(regions));
  regions.forEach ( region => {
    region.start = region.start || 0;
    // region.end is being looked up and filled in mapper
    // and thus it should exist by this point.
    assert(region.end);
  });
  return;
};

let reference2RegionString = (queryReference, refSepEncoding) => {
  assert(typeof refSepEncoding !== 'undefined');
  let refSeparator;
  let regionString = queryReference.referenceName;
  var _getInt = (str) => {
    /* Suppress unhelpful advice not to use Number as a constructor */
    /* jshint -W053 */
    let n = new Number(str);
    /* jshint +W053 */
    n = n.valueOf();
    if (!Number.isInteger(n)) {
      throw new RangeError('ranges must be integers');
    }
    if (n < 0) {
      throw new RangeError('ranges must be unsigned integers');
    }
    return n;
  };
  assert(regionString.match(/^[-\w\.]+$/),
         'Invalid character in reference name ' + queryReference.referenceName);
  if (refSepEncoding) {
    refSeparator = encodeURIComponent(':');
  } else {
    refSeparator = ':';
  }
  /*
   * GA4GH ranges are semi-open starting from zero.
   * Samtools ranges are closed starting from one.
   */
  let start = queryReference.start || '0';
  start = _getInt(start);
  let end = (typeof queryReference.end === 'undefined') ? null : _getInt(queryReference.end);
  if (end !== null && end <= start) {
    throw new RangeError('Range end should be bigger than start');
  }
  start += 1;

  if (end) {
    regionString += (refSeparator + start + '-' + end);
  } else {
    if (start !== 1) {
      regionString += refSeparator + start;
    }
  }
  return regionString;
};

/** Application model */
class RangerModel {

  /**
   * Creates a RangerModel object instance.
   */
  constructor() {}

  /**
   * Retrieves, processes and streams the data.
   *
   * For a single file query, retrieves one data source (file),
   * optionally selects the requested regions, converts the
   * output to the requested format and streams data to destination.
   *
   * For a miltiple files query, retrieves multiple data
   * source(files), merges them, optionally selects the requested region,
   * marks duplicates and streams data to destination.
   *
   * @param query       - an object representing a request
   * @param destination - a writable stream for the output
   * @param endResponseClback - callback to execute at the end streaming
   */
  process(query, destination, endResponseClback) {
    assert(query, 'Query object is required');
    assert.equal(typeof query, 'object', 'Query should be an object');
    assert.equal(typeof query.files, 'object', "Query should have a 'files' object inside");
    assert.equal(typeof query.files.length, 'number', "Files object should have numeric 'length' property");
    assert.notEqual(query.files.length, 0, 'Files should be given');
    assert(destination, 'Destination stream is required');
    assert ((typeof endResponseClback == 'function'), 'End callback is required');
    if (query.files.length === 1) {
      this._getFile(query, destination, endResponseClback);
    } else {
      this._mergeFiles(query, destination, endResponseClback);
    }
  }

  /**
   * Returns default data format - bam. Class method.
   */
  static defaultFormat() {
    return DEFAULT_FORMAT;
  }

  /**
   * Returns an array containing supported data formats. Class method.
   */
  static supportedFormats() {
    return [DEFAULT_FORMAT, 'CRAM', 'SAM', 'VCF'];
  }

  /**
   * Returns a boolean result indicating whether teh given format
   * is supported. Class method.
   */
  static supportsFormat(format) {
    assert(format, 'Non-empty format string should be given');
    return this.supportedFormats().indexOf(format) >= 0;
  }

  /**
   * Returns an array containing readable text data formats. Class method.
   */
  static textualFormats() {
    return ['SAM', 'VCF'];
  }

  /**
   * Returns an array containing readable text data formats. Class method.
   */
  static isTextualFormat(format) {
    assert(format, 'Non-empty format string should be given');
    return this.textualFormats().indexOf(format) >= 0;
  }

  /**
   * Returns true if the query has a reference property. Class method.
   */
  static hasReference(query) {
    assert(query, 'Query must be given');
    return !!query.reference;
  }

  _stViewAttrs(query) {
    assert(query);
    let attrs = ['view', '-h'];
    let format = query.format;
    if ( format ) {
      // To create .vcf, Freebayes needs input to be in .bam format.
      // So we will first convert cram from database into .bam, then
      // process in Freebayes.
      if ( format === 'BAM' || format === 'VCF' ) {
        attrs.push('-b');
      } else if ( format === 'CRAM' ) {
        attrs.push('-C');
      }
      attrs.push('-T', RangerModel.fixReference(query));
    }
    attrs.push(query.files.shift() || "-");
    if (query.region) {
      attrs = attrs.concat(query.region);
    } else if (query.multiregions_key) {
      query.multiregions_key.forEach( region => {
        attrs = attrs.concat(reference2RegionString(region, false));
      });
    }
    LOGGER.debug('view attrs: ' + attrs);
    return attrs;
  }

  _stMergeAttrs(query, dir) {
    assert(query);
    let attrs = ['merge', "-u"];
    if (query.region) {
      attrs.push('-R', query.region);
    } else if (query.multiregions_key) {
      // Changes the value of query's region.start due to JS shenanigans.
      checkStartEnd(query.multiregions_key);
      let bedString = buildBedFile(query.multiregions_key);
      let fileName = Date.now() + BED_FILETYPE;
      fse.ensureDirSync(dir + TEMPDIR_FILEPATH);
      fse.writeFileSync(dir + TEMPDIR_FILEPATH + '/' + fileName, bedString);
      attrs.push('-L');
      attrs.push(dir + TEMPDIR_FILEPATH + '/' + fileName);
    }
    attrs.push('--reference', RangerModel.fixReference(query));
    attrs.push('-');
    let files = query.files;
    let re_bam    = /\.bam$/i;  // Ignore case with /i
    let re_cram   = /\.cram$/i; // Ignore case with /i
    let some_bam  = files.some(function(f) { return re_bam.test(f.data_object); });
    let some_cram = files.some(function(f) { return re_cram.test(f.data_object); });
    if (some_bam && some_cram) {
      throw new Error('Inconsistent format, all files should be either bam or cram');
    }
    attrs = attrs.concat(files);
    LOGGER.debug('merge attrs: ' + attrs);
    return attrs;
  }

  static fixReference(query) {
    assert(query);
    let options = config.provide();
    let rootref = options.get('references');
    rootref = rootref ? rootref + '/' : '';
    return path.resolve(rootref + query.reference);
  }

  _fbAttrs(query) {
    assert(query);
    assert(RangerModel.hasReference(query),
      'database does not hold location of reference .fa file');
    // '-c' option allows freebayes to read from stdin ie. from pipe
    let attrs = ['-c'];
    attrs.push('-f', RangerModel.fixReference(query));
    if (query.region) {
      attrs.push('-r', query.region);
    }
    LOGGER.debug('freebayes attrs: ' + attrs);
    return attrs;
  }

  _scheduleVCF(query, processes, dir) {
    const freebayes = child.spawn(
      FREEBAYES_COMMAND,
      this._fbAttrs(query),
      {cwd: dir}
    );
    freebayes.title = 'freebayes';
    processes.push(freebayes);
  }

  _getFile(query, destination, endResponseClback) {
    let dir = config.tempFilePath();
    fse.ensureDirSync(dir);

    const view = child.spawn(
      SAMTOOLS_COMMAND,
      this._stViewAttrs(query),
      { cwd: dir }
    );
    view.title = 'samtools view';
    let processes = [view];
    if (query.format === 'VCF') {
      this._scheduleVCF(query, processes, dir);
    }
    pipeline(
      processes,
      () => { endResponseClback(false); _cleanup( dir ); },
      () => { endResponseClback(true);  _cleanup( dir ); } )
    .run(destination);
  }

  _mergeFiles(query, destination, endResponseClback) {
    let dir = config.tempFilePath();
    fse.ensureDirSync(dir);
    let attrs = this._stMergeAttrs(query, dir);
    const merge = child.spawn(
      SAMTOOLS_COMMAND,
      attrs,
      {cwd: dir}
    );
    merge.title = 'samtools merge';

    /* The full, heavy Samtools Pipeline is below,
       which is disabled by default and not supported by us.
       To enable this, uncomment this block and another marked
       section below.
       We use the light Pipeline as the full version may be too
       computationally taxing, and the majority of our own data can
       work with the light version. (Older data is the exception.)
       Explained in the README.

    // -O to output to stdout, - to get data in from stdin
    attrs = ['collate', '-O', '-u', '-'];
    const collate = child.spawn(
      SAMTOOLS_COMMAND,
      attrs,
      {cwd: dir}
    );
    collate.title = 'samtools collate';

    // two '-', one to pipe in from stdin and one to pipe to stdout
    // -m to add ms tags for markdup
    attrs = ['fixmate', '-m', '-', '-'];
    const fixmate = child.spawn(
      SAMTOOLS_COMMAND,
      attrs,
      {cwd: dir}
    );
    fixmate.title = 'samtools fixmate';

    // pipe output to stdout automatically, - to pipe from stdin
    attrs = ['sort', '-'];
    const sort = child.spawn(
      SAMTOOLS_COMMAND,
      attrs,
      {cwd: dir}
    );
    sort.title = 'samtools sort';
    */

    // two '-', one to pipe in from stdin and one to pipe to stdout
    attrs = ['markdup', '-', '-'];
    const markdup = child.spawn(
      SAMTOOLS_COMMAND,
      attrs,
      {cwd: dir}
    );
    markdup.title = 'samtools markdup';

    // multiregions_key is removed as it would have already been used
    // during the merge above. If it is passed through into the view
    // post-merge, it'll fail the process as it is trying to operate on
    // an unindexed file.
    delete query.multiregions_key;
    delete query.region;
    delete query.directory;
    query.files = [];
    const view  = child.spawn(
      SAMTOOLS_COMMAND,
      this._stViewAttrs(query),
      {cwd: dir}
    );
    view.title = 'samtools view (post-merge)';
    let processes = [merge, markdup, view];
    // Below is the full Samtools Pipeline. Disabled as the light version
    // is used by default. More information is above, and in the README.
    // let processes = [merge, collate, fixmate, sort, markdup, view];
    if (query.format === 'VCF') {
      this._scheduleVCF(query, processes, dir);
    }

    pipeline(
      processes,
      () => { endResponseClback(false); _cleanup( dir ); },
      () => { endResponseClback(true);  _cleanup( dir ); } )
    .run(destination);
  }
}

module.exports.RangerModel = RangerModel;
module.exports.reference2RegionString = reference2RegionString;
module.exports.buildBedFile = buildBedFile;
module.exports.checkStartEnd = checkStartEnd;
