"use strict";

const LOGGER = require('winston');
const moment = require('moment');

var getTimeStamp = function() {
  return moment().format('YYYY-MM-DD HH:mm:ss Z');
};

try {
  LOGGER.remove(LOGGER.transports.Console);
  LOGGER.add(LOGGER.transports.Console, {
    level:     'debug',
    colorize:  true,
    json:      false,
    timestamp: getTimeStamp,
    stderrLevels: [ 'error' ] // Which levels to send to stderr
  });
  LOGGER.add(LOGGER.transports.File, {
    level:    'debug',
    filename: 'npg_ranger.winston.log',
    timestamp: getTimeStamp,
    maxsize:   10 * 1000 * 1000,
    maxFiles:  20
  });
} catch (e) {
  console.error('Error when trying to init logging with Winston');
  console.error(e);
  process.exit(1);
}

module.exports = LOGGER;
