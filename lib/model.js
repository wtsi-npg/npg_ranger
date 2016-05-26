"use strict";

/**
 * model module.
 * @module lib/model
 *
 * @description Model for the streaming HTTP server.
 *
 * @example <caption>Example usage of the model.</caption>
 *   const RangerModel = require('../lib/model.js');
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
const os       = require('os');
const path     = require('path');
const pipeline = require('../lib/pipeline.js');

const SAMTOOLS_COMMAND     = 'samtools';
const BBB_MARKDUPS_COMMAND = 'bamstreamingmarkduplicates';
const DEFAULT_FORMAT       = 'bam';

/** Application model */
class RangerModel {

  /**
   * Creates a RangerModel object instance.
   * @param tmpDir - an optional path to a directory for temporary data
   *                 defaults to the OS tmp directory
   */
  constructor(tmpDir) {
    this.tmpDir = (tmpDir || os.tmpdir());
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
    assert(destination, 'Destination stream is required');
    assert ((typeof endResponseClback == 'function'), 'End callback is required');
    if (query.files.length === 1) {
      this._getFile(query, destination, endResponseClback);
    } else {
      this._mergeFiles(query, destination, endResponseClback);
    }
  }

  /**
   * Returns default data format - bam.
   */
  defaultFormat() {
    return DEFAULT_FORMAT;
  }

  _tempFilePath() {
    return path.join(this.tmpDir, Math.random().toString().substr(2));
  }

  _stViewAttrs(query) {
    assert(query);
    let attrs = ['view', '-h'];
    if (query.format && (query.format === 'bam' || query.format === 'cram')) {
      attrs.push(query.format === 'bam' ? '-b' : '-C');
    }
    attrs.push(query.files.shift() || "-");
    if (query.region) {
      attrs = attrs.concat(query.region);
    }
    console.log('view attrs: ' + attrs);
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
    let re_bam    = /\.bam$/;
    let re_cram   = /\.cram$/;
    let some_bam  = files.some(function(f) { return re_bam.test(f.data_object); });
    let some_cram = files.some(function(f) { return re_cram.test(f.data_object); });
    if (some_bam && some_cram) {
      throw new Error('Inconsistent format, all files should be either bam or cram');
    }

    attrs = attrs.concat(files);
    console.log('merge attrs: ' + attrs);
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
          console.log(`Failed to remove ${dir}: ${err}`);
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