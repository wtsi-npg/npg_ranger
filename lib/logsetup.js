"use strict";

const winston = require('winston');
const moment = require('moment');

/**
 * @external winston
 * @see      {@link https://www.npmjs.com/package/winston|winston}
 */

/**
 * @external moment
 * @see      {@link https://www.npmjs.com/package/moment|moment}
 */

/**
 * @module logsetup
 *
 * @requires {@link external:winston|winston}
 * @requires {@link external:moment|moment}
 *
 * @copyright Genome Research Limited 2022
 */

var getTimeStamp = function() {
  return moment().format('YYYY-MM-DD HH:mm:ss Z');
};


const LOGGER = new winston.createLogger({
  transports: [
    new winston.transports.Console( {
      colorize:  true,
      json:      false,
      timestamp: getTimeStamp,
      stderrLevels: [ 'error', 'warn', 'info', 'verbose', 'debug', 'silly' ] // Which levels to send to stderr (all)
    })
  ]
});

const BINLOGGER = new winston.createLogger({
  level: 'error',
  transports: [
    new winston.transports.Console( {
      level:     'error',
      colorize:  true,
      json:      false,
      timestamp: getTimeStamp,
      stderrLevels: [ 'error', 'warn', 'info', 'verbose', 'debug', 'silly' ] // Which levels to send to stderr (all)
    })
  ]
});


module.exports = {LOGGER, BINLOGGER};
