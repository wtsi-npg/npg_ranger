"use strict";

/**
 * mapper module.
 * @module lib/mapper
 *
 * @description Authorisation for access to data sources (sequencing files).
 *              Relies on data access group being known for each data
 *              source and the user having access to a whole set or subset
 *              of data access groups.
 *
 * @example <caption>Example usage of the mapper module.</caption>
 *   const DataAccess = require('../lib/mapper.js');
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
 *   da.on('data', (files) => {
 *          console.log('Data ready');
 *     // Do something with data
 *     // This is how files array might look
 *     // [ {"filepath_by_host" : {"*" : "irods:/seq/foo1"},
 *     //    "access_control_group_id": "6"},
 *     //   {"filepath_by_host" : {"*" : "irods:/seq/foo2"},
 *     //    "access_control_group_id": "7"}]
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

const EventEmitter = require('events');

const FILE_INFO_COLLECTION  = 'fileinfo';
const ERROR_EVENT_NAME      = 'error';
const NO_DATA_EVENT_NAME    = 'nodata';
const DATA_READY_EVENT_NAME = 'data';

/** Class mapping client query to the location and access info of sequencing files. */
class DataMapper extends EventEmitter {
  /**
   * Creates an DataMapper type object.
   * @param db - mpngodb connection object
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
   * Retirns an array of objects with information about sequencing files.
   * @param query - an object having either name or accession attribute defined
   * @param host  - host name
   */
  getFileInfo(query, host) {
    if (!query) {
      throw new ReferenceError('Query hash is required');
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
    console.log(dbquery);

    let columns = { _id: 0,
                    access_control_group_id: 1,
                    'filepath_by_host.*': 1 };
    columns['filepath_by_host.' + host] = 1;

    var cursor = this.db.collection(FILE_INFO_COLLECTION).find(dbquery, columns);
    var files   = [];
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
            this.emit(NO_DATA_EVENT_NAME,
              'No files for ' + (a ? a : query.name));
          } else {
            let data = {};
            data.files        = files.map((f) => {
              return f.filepath_by_host[host] || f.filepath_by_host["*"];
            });
            data.accessGroups = files.map((f) => { return f.access_control_group_id; });
            this.emit(DATA_READY_EVENT_NAME, data);
          }
        }
      }
    });
  }
}

module.exports = DataMapper;

