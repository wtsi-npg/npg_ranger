"use strict";

const LOGGER = require('winston');
const moment = require('moment');

const config = require('../lib/config.js');

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
 * @author Jaime Tovar
 * @copyright Genome Research Limited 2016
 */

var getTimeStamp = function() {
  return moment().format('YYYY-MM-DD HH:mm:ss Z');
};

var setup = function() {
  try {
    LOGGER.remove(LOGGER.transports.Console);
    if (config.provide().get('logconsole')) {
      LOGGER.add(LOGGER.transports.Console, {
        level:     'warn',
        colorize:  true,
        json:      false,
        timestamp: getTimeStamp,
        stderrLevels: [ 'error' ] // Which levels to send to stderr
      });
    }
    LOGGER.add(LOGGER.transports.File, {
      level:    'warn',
      filename: 'npg_ranger.winston.log',
      json:      false,
      timestamp: getTimeStamp,
      maxsize:   10 * 1000 * 1000,
      maxFiles:  20
    });
  } catch (e) {
    console.error('Error when trying to init logging with Winston');
    console.error(e);
    process.exit(1);
  }
};

module.exports = LOGGER;
module.exports.setup = setup;
