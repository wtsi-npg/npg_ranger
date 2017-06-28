"use strict";

const TOKEN_BEARER_KEY_NAME = 'Authorization';
const TOKEN_CONFIG_KEY_NAME = 'token';

const AUTH_TYPE_TOKEN = 'authtoken';
const AUTH_TYPE_USER  = 'authuser';
const AUTH_URL_TOKEN  = '/validateToken';
const AUTH_URL_USER   = '/validateUser';

const DEFAULT_PROTOCOL = 'http:';
const SECURE_PROTOCOL  = 'https:';

module.exports = {
  TOKEN_BEARER_KEY_NAME,
  TOKEN_CONFIG_KEY_NAME,
  AUTH_TYPE_TOKEN,
  AUTH_TYPE_USER,
  AUTH_URL_TOKEN,
  AUTH_URL_USER,
  DEFAULT_PROTOCOL,
  SECURE_PROTOCOL
};
