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
 * @copyright Genome Research Limited 2017
 */

var getTimeStamp = function() {
  return moment().format('YYYY-MM-DD HH:mm:ss Z');
};

const LOGGER = winston.createLogger({
  level:     'error',
  transports: [
    new winston.transports.Console()
  ]
});

LOGGER.configure({
  colorize:  true,
  json:      false,
  stderrLevels: [ 'error', 'warn', 'info', 'verbose', 'debug', 'silly' ], // Which levels to send to stderr (all)
  timestamp: getTimeStamp,
});

module.exports = LOGGER;
