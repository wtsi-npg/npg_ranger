"use strict";

/**
 * auth module.
 * @module lib/auth
 *
 * @description Authorisation for access to data sources (sequencing files).
 *              Relies on data access group being known for each data
 *              source and the user having access to a whole set or subset
 *              of data access groups.
 *
 * @example <caption>Example usage of the auth module.</caption>
 *   const DataAccess = require('../lib/auth.js');
 *   // Define/get info about data sources.
 *   var files = [
 *     {"filepath_by_host" : {"*" : "irods:/seq/foo1"},
 *         "access_control_group_id": "6"},
 *    {"filepath_by_host" : {"*" : "irods:/seq/foo2"},
 *         "access_control_group_id": "7"}
 *               ];
 *
 *   // Example of a database record mapping users to access groups
 *   var records_json = '[
 *     {
 *       "access_control_group_id": "6",
 *       "members": ["alice", "bob", "charles"]
 *     },
 *     {
 *       "access_control_group_id": "7",
 *       "members": ["alice", "bob", "james"]
 *     }
 *   ]';
 *   // Create a new object.
 *   var da = new DataAccess(db, files);
 *   // Register listeners for events this object emits.
 *   da.on('authorised', (username) => {
 *          console.log(`User ${username} is given access`);
 *        });
 *   da.on('failed', (username, message) => {
 *          console.log(
 *            `Authorisation failed for user '${username}': ${message}`);
 *        });
 *   // Call authorisation for a user, giving a username
 *   da.authorise('alice');
 *
 * @author Marina Gourtovaia
 * @copyright Genome Research Limited 2016
 */


const EventEmitter = require('events');

const ACCESS_CONTROL_COLLECTION = 'access_control_group';
const AUTH_EVENT_NAME           = 'authorised';
const AUTH_FAILED_EVENT_NAME    = 'failed';

/** Class encupsulating authorisation. */
class DataAccess extends EventEmitter {
  /**
   * Creates an DataClass type object.
   * @param db - mpngodb connection object
   * @param files - an array of hashes describing a data source,
   *                see example for more details
   */
  constructor(db, files) {
    super();
    if (!db) {
      throw new ReferenceError('Database handle is required');
    }
    if (!files || (files instanceof Array === false) || (files.length === 0)) {
      throw new Error('File info is not available');
    }
    this.db    = db;
    this.files = files;
  }

  /**
   * Authorise the user and emit an appropriate event.
   * @param username - username as a string
   */
  authorise(username) {
    if (!username) {
      throw new ReferenceError('Username is required');
    }
    var agroup_ids = this.files.map(function(file) { return file.access_control_group_id; });
    if (!agroup_ids.every(function(id) { return id; })) {
      this.emit(AUTH_FAILED_EVENT_NAME, username,
        'Some of files do not have access group defined');
      return;
    }
    agroup_ids.sort();
    agroup_ids = agroup_ids.filter(function(item, index, thisArray) {
      return (index === 0) ? 1 : ((item === thisArray[index - 1]) ? 0 : 1);
    });
    console.log('ACCESS GROUP IDS: ' + agroup_ids.join(' '));
    var dbquery = {
      members:                 username,
      access_control_group_id: {$in: agroup_ids}
    };

    var count_callback = (err, count) => {
      if (err) {
        this.emit(AUTH_FAILED_EVENT_NAME, username,
          'Failed to get authorisation info, DB error: ' + err);
      } else {
        if (count == agroup_ids.length) {
          this.emit(AUTH_EVENT_NAME, username);
        } else {
          this.emit(AUTH_FAILED_EVENT_NAME, username,
            'Not authorised for ' + (count ? 'some' : 'any') + ' of the files');
        }
      }
    };
    this.db.collection(ACCESS_CONTROL_COLLECTION).count(dbquery, count_callback);
  }
}

module.exports = DataAccess;

