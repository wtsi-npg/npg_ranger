"use strict";

const fs = require('fs');

/**
 * @module test/server/test_utils
 *
 * @description
 * <p>Common functions for test modules</p>
 *
 * @requires {@link external:fs|fs}
 *
 * @author Jaime Tovar
 * @copyright Genome Research Limited 2017
 */

/**
 * Best-effort approach to remove socket files from the system. If something
 * goes wrong it will report the error to console.
 * @param  {string} socket socket path
 */
var removeSocket = ( socket ) => {
  try {
    // check if constants is defined for compatibility between node v4 and v6.3+
    let flag = fs.constants ? fs.constants.W_OK : fs.W_OK;
    fs.access(socket, flag, ( err ) => {
      if ( !err ) {
        fs.unlinkSync(socket);
      } else {
        console.log( err );
      }
    });
  } catch (e) { console.log(e); }
};

module.exports = {
  removeSocket: removeSocket
};
