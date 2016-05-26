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
 *   let m = new RangerModel(
 *      response, query, errorResponseClback, endResponseClback);
 *   // Create a new object and define the temporary data directory.
 *   let m1 = new RangerModel(
 *      response, query, errorResponseClback, endResponseClback, dir);
 *   // stream data for one file
 *   m.getFile();
 *   // stream merged data for multiple files
 *   m.mergeFiles();
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

const SAMTOOLS_COMMAND        = 'samtools';
const BBB_MARKDUPS_COMMAND    = 'bamstreamingmarkduplicates';

/** Application model */
class RangerModel {

  /**
   * Creates a RangerModel object instance.
   * @param response - HTTP response object
   * @param query    - an object representing a request
   * @param errorResponseClback -
   *   callback to invoke to create an HTTP error response
   * @param endResponseClback -
   *   callback to invoke to close the response
   * @param tmpDir   - an optional path to a directory for
   *                   for temporary data, defaults to the
   *                   OS tmp directory
   */
  constructor(
    response, query, errorResponseClback, endResponseClback, tmpDir) {

    assert(response, 'HTTP response object is required');
    this.response = response;

    assert(query, 'Query object is required');
    this.query = query;

    assert((typeof errorResponseClback == 'function'),
      'Error callback is required');
    this.errorResponseClback = errorResponseClback;
    assert ((typeof endResponseClback == 'function'),
      'End callback is required');
    this.endResponseClback = endResponseClback;

    this.tmpDir = (tmpDir || os.tmpdir());
  }

  tempFilePath() {
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
    attrs.push('tmpfile=' + this.tempFilePath());
    attrs.push('M=/dev/null');
    return attrs;
  }

  /**
   * Retrieves one data source(file), optionally selects
   * the requested regions, converts the output to the
   * requested format and streams data to destination.
   */
  getFile() {
    const view = child.spawn(SAMTOOLS_COMMAND, this._stViewAttrs(this.query));
    view.title = 'samtools view';
    pipeline(
      [view],
      () => {this.endResponseClback(this.response,false);},
      () => {this.endResponseClback(this.response,true);} )
    .run(this.response);
  }

  /**
   * Retrieves multiple data source(files), merges them,
   * optionally selects the requested region, marks
   * duplicates and streams data to destination.
   */
  mergeFiles() {

    let dir = this.tempFilePath();
    fs.mkdirSync(dir);

    const cleanup = function() {
      fse.remove(dir, function(err) {
        if (err) {
          console.log(`Failed to remove ${dir}: ${err}`);
        }
      });
    };

    var attrs;
    try {
      attrs = this._stMergeAttrs(this.query);
    } catch (ex) {
      this.errorResponseClback(this.response, 500, ex);
    }
    if (!attrs) {
      return;
    }

    const merge = child.spawn(SAMTOOLS_COMMAND,
                  attrs,
                  {cwd: dir});
    merge.title = 'samtools merge';

    const markdup = child.spawn(BBB_MARKDUPS_COMMAND, this._bbbMarkDupsAttrs());
    markdup.title = BBB_MARKDUPS_COMMAND;

    delete this.query.region;
    delete this.query.directory;
    const view  = child.spawn(SAMTOOLS_COMMAND, this._stViewAttrs(this.query));
    view.title = 'samtools view (post-merge)';

    pipeline(
      [merge,markdup,view],
      () => {this.endResponseClback(this.response,false); cleanup();},
      () => {this.endResponseClback(this.response,true);  cleanup();} )
    .run(this.response);
  }
}

module.exports = RangerModel;