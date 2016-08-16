"use strict";

/**
 * mapper module.
 * @module server/mapper
 *
 * @description Authorisation for access to data sources (sequencing files).
 *              Relies on data access group being known for each data
 *              source and the user having access to a whole set or subset
 *              of data access groups.
 *
 * @example <caption>Example usage of the mapper module.</caption>
 *   const DataAccess = require('../lib/server/mapper.js');
 *   // Create a new object.
 *   var da = new DataMapper(db);
 *   // Register listeners for events this object emits.
 *   da.on('error', (err) => {
 *          console.log(`Error retrieving data: ${err}`);
 *   });
 *   da.on('nodata', (message) => {
 *          console.log(
 *            `No data retrieved: ${message}`);
 *   });
 *   da.on('data', (data) => {
 *          console.log('Data ready');
 *     // Do something with data
 *     // This is how the data object might look
 *     // [{file: "irods:/seq/foo1", accessGroups: "6"},
 *     //  {file: "irods:/seq/foo2", accessGroups: "7"}]
 *   });
 *   // Call mapper
 *   let hostname = 'localhost';
 *   let query = {'accession': 'XYZ45678'};
 *   // Other options:
 *   // let query = {'name': 'my.cram'};
 *   // let query = {'name': 'my.cram', 'directory': 'some'};
 *   dm.getFileInfo(query, hostname);
 *
 * @author Marina Gourtovaia
 * @copyright Genome Research Limited 2016
 */

const LOGGER       = require('winston');

const EventEmitter = require('events');

const FILE_INFO_COLLECTION  = 'fileinfo';
const ERROR_EVENT_NAME      = 'error';
const NO_DATA_EVENT_NAME    = 'nodata';
const DATA_READY_EVENT_NAME = 'data';

/** Class mapping client query to the location and access info of sequencing files. */
class DataMapper extends EventEmitter {
  /**
   * Creates a DataMapper type object.
   * @param db - mongodb connection object
   */
  constructor(db) {
    super();
    if (!db) {
      throw new ReferenceError('Database handle is required');
    }
    this.db = db;
  }

  /**
   * Maps the query to the location and access info of sequencing files.
   * Tries to find files co-located with the given host.
   * Returns an array of objects with information about sequencing files.
   * @param query - an object having either name or accession attribute defined
   * @param host  - host name
   */
  getFileInfo(query, host) {
    if (!query) {
      throw new ReferenceError('Query object is required');
    }
    if (!host) {
      throw new ReferenceError('Host name is required');
    }

    var a = query.accession;
    var dbquery;
    if (a) {
      dbquery =  { $and: [{'avh.sample_accession_number': a},
                          {'avh.target':    "1"},
                          {'avh.manual_qc': "1"},
                          {'avh.alignment': "1"} ]};
    } else if (query.name) {
      dbquery =  {data_object: query.name};
      if (query.directory) {
        dbquery =  { $and: [ dbquery, {collection: query.directory} ]};
      }
    } else {
      throw new Error('Sample accession number or file name should be given');
    }

    LOGGER.debug(dbquery);

    let columns = { _id: 0,
                    access_control_group_id: 1,
                    'filepath_by_host.*': 1,
                    'avh.reference': 1 };
    columns['filepath_by_host.' + host] = 1;

    var cursor = this.db.collection(FILE_INFO_COLLECTION).find(dbquery, columns);
    var files   = [];
    let noFilesMessage = 'No files for ' + (a ?
      `sample accession ${a}` :
      query.name + (query.directory ? ' in ' + query.directory : '')
    );
    cursor.each( (err, doc) => {
      if (err) {
        try {
          cursor.close();
        } catch (e) {}
        this.emit(ERROR_EVENT_NAME,
          'Failed to map input to files, DB error: ' + err);
      } else {
        if (doc != null) {
          files.push(doc);
        } else {
          cursor.close();
          if (files.length === 0) {
            this.emit(NO_DATA_EVENT_NAME, noFilesMessage);
          } else {
            let data = files.map( (f) => {
              LOGGER.info(f);
              let d = {};
              d.file        = f.filepath_by_host[host]  || f.filepath_by_host["*"];
              d.accessGroup = f.access_control_group_id || '';

              if (f.avh.reference) {
                // Need to change the reference path to ensure correctness as
                // database may be incorrect. Assume the last 4 folder names
                // are correct (except the very last).
                let ref   = f.avh.reference.split('/').slice(-5);
                ref[ref.length - 2] = 'fasta';
                ref                 = ref.join('/');
                d.reference         = '/lustre/scratch110/srpipe/references/' + ref;
                LOGGER.info('mapper thinks refpath looks like: ');
                LOGGER.info(d.reference);
              }

              return d;
            });
            data = data.filter( (d) => {return d.file;} );
            if (data.length) {
              this.emit(DATA_READY_EVENT_NAME, data);
            } else {
              this.emit(NO_DATA_EVENT_NAME, noFilesMessage);
            }
          }
        }
      }
    });
  }
}

module.exports = DataMapper;

