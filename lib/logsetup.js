"use strict";

const LOGGER = require('winston');

try {
  LOGGER.remove(LOGGER.transports.Console);
  LOGGER.add(LOGGER.transports.Console, {
    level:        'debug',
    colorize:     true,
    timestamp:    true,
    json:         false,
    stderrLevels: [ 'error' ]
  });
  LOGGER.add(LOGGER.transports.File, {
    level:    'debug',
    filename: 'npg_ranger.winston.log',
    maxsize:  10 * 1000 * 1000,
    maxFiles: 20
  });
} catch (e) {
  console.log('Error when trying to init winston');
  console.log(e);
  process.exit(1);
}

module.exports = LOGGER;
