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
 *   // Create a new object and define the temporary data directory.
 *   m = new RangerModel(dir);
 *   // process data
 *   m.process(query, destination, endResponseClback);
 *
 * @author Marina Gourtovaia
 * @copyright Genome Research Limited 2016
 */

const assert   = require('assert');
const child    = require('child_process');
const fs       = require('fs');
const fse      = require('fs-extra');
const path     = require('path');
const LOGGER   = require('winston');
const pipeline = require('./pipeline.js');
const config   = require('../config.js');

const SAMTOOLS_COMMAND     = 'samtools';
const BBB_MARKDUPS_COMMAND = 'bamstreamingmarkduplicates';
const DEFAULT_FORMAT       = 'BAM';

/** Application model */
class RangerModel {

  /**
   * Creates a RangerModel object instance.
   * @param tmpDir - an optional path to a directory for temporary data
   *                 defaults to the OS tmp directory
   */
  constructor(tmpDir) {
    var options = config.provide();
    this.tmpDir = tmpDir || options.get('tempdir');
    assert(fs.existsSync(this.tmpDir),
      `Temp data directory '${this.tmpDir}' does not exist`);
  }

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
    return [DEFAULT_FORMAT, 'CRAM', 'SAM'];
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
    return ['SAM'];
  }

  /**
   * Returns an array containing readable text data formats. Class method.
   */
  static isTextualFormat(format) {
    assert(format, 'Non-empty format string should be given');
    return this.textualFormats().indexOf(format) >= 0;
  }

  _tempFilePath() {
    return path.join(this.tmpDir, Math.random().toString().substr(2));
  }

  _stViewAttrs(query) {
    assert(query);
    let attrs = ['view', '-h'];
    let format = query.format;
    if ( format ) {
      if ( format === 'BAM' || format === 'CRAM' ) {
        attrs.push(format === 'BAM' ? '-b' : '-C');
      }
    }
    attrs.push(query.files.shift() || "-");
    if (query.region) {
      attrs = attrs.concat(query.region);
    }
    LOGGER.debug('view attrs: ' + attrs);
    return attrs;
  }

  _stMergeAttrs(query) {
    assert(query);
    let attrs = ['merge', "-u"];
    if (query.region) {
      let regions = query.region;
      if (typeof regions != 'object') {
        regions = [regions];
      }
      regions.map(function(r) {
        attrs.push('-R');
        attrs.push(r);
      });
    }
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

  _bbbMarkDupsAttrs() {
    var attrs = ['level=0','verbose=0','resetdupflag=1'];
    attrs.push('tmpfile=' + this._tempFilePath());
    attrs.push('M=/dev/null');
    return attrs;
  }

  _getFile(query, destination, endResponseClback) {
    const view = child.spawn(SAMTOOLS_COMMAND, this._stViewAttrs(query));
    view.title = 'samtools view';
    pipeline(
      [view],
      () => {endResponseClback(false);},
      () => {endResponseClback(true);} )
    .run(destination);
  }

  _mergeFiles(query, destination, endResponseClback) {
    let dir = this._tempFilePath();
    fs.mkdirSync(dir);

    const cleanup = function() {
      fse.remove(dir, function(err) {
        if (err) {
          LOGGER.warn(`Failed to remove ${dir}: ${err}`);
        }
      });
    };

    let attrs = this._stMergeAttrs(query);
    const merge = child.spawn(SAMTOOLS_COMMAND,
                  attrs,
                  {cwd: dir});
    merge.title = 'samtools merge';

    const markdup = child.spawn(BBB_MARKDUPS_COMMAND, this._bbbMarkDupsAttrs());
    markdup.title = BBB_MARKDUPS_COMMAND;

    delete query.region;
    delete query.directory;
    query.files = [];
    const view  = child.spawn(SAMTOOLS_COMMAND, this._stViewAttrs(query));
    view.title = 'samtools view (post-merge)';

    pipeline(
      [merge,markdup,view],
      () => {endResponseClback(false); cleanup();},
      () => {endResponseClback(true);  cleanup();} )
    .run(destination);
  }
}

module.exports = RangerModel;
