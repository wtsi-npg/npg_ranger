"use strict";

/**
 * Helper module to stablish single point of entry to the package
 * @module main
 *
 * @author Jaime Tovar
 * @copyright Genome Research Limited 2016
 */

/**
 * Business logic abstraction of GA4GH and npg_ranger requests
 * @type {RangerRequest}
 */
exports.RangerRequest = require('./client/rangerRequest').RangerRequest;
exports.uiclient      = require('./uiclient/accUtils');
