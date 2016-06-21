"use strict";

// For http
const EventEmitter = require('events');

// const assert       = require('assert');

// const http         = require('http');
// const url          = require('url');

class HTTPTransport extends EventEmitter {

  constructor () {
    super();
  }

  doNothing() {
    return;
  }
}

module.exports = HTTPTransport;
