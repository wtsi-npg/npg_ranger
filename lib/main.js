"use strict";

/**
 * Helper module to stablish single point of entry to the package
 * @module main
 *
 * @copyright Genome Research Limited 2017
 */

/**
 * Business logic abstraction of GA4GH and npg_ranger requests
 * @type {RangerRequest}
 */
exports.RangerRequest = require('./client/rangerRequest').RangerRequest;
exports.uiclient      = require('./uiclient/accUtils');
