"use strict";

const fs = require('fs');

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
