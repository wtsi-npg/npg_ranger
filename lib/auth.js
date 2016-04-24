"use strict";

const EventEmitter = require('events');

const ACCESS_CONTROL_COLLECTION = 'access_control_group';
const AUTH_EVENT_NAME           = 'authorised';
const AUTH_FAILED_EVENT_NAME    = 'failed';

class DataAccess extends EventEmitter {

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
      console.log('COUNT ' + count);
      console.log('ERROR ' + err);
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

