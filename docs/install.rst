=====================
Installing npg_ranger
=====================

`When cloning from github`_


`Install with npm`_

When cloning from github
========================

Clone the project
-----------------

::

  git clone -b master https://github.com/wtsi-npg/npg_ranger.git && cd npg_ranger

Building and installing
=======================

Install node dependencies
-------------------------

::

  npm install

Install with npm
================

For local installation
----------------------

Will install npg_ranger locally in the current folder. Links for executables (npg_ranger_server and npg_ranger_client) can be found in ./node_modules/.bin or directly the executables (server.js and client.js) as part of the package structure in ./node_modules/npg_ranger/bin

::

  npm install npg_ranger

For global installation
-----------------------

This will install npg_ranger globally for your node deployment. Unless you installed nodejs locally is likely you will need sudo to use this option. It will modify your path to include aliases for npg_ranger executables (npg_ranger_server and npg_ranger_client).

::

  npm install -g npg_ranger

Post install
============

Generate documentation
----------------------

::

  grunt jsdoc

Test
----

::

  grunt lint                           - runs lint test
  grunt jasmine                        - runs jasmine tests
  grunt jasmine -v                     - runs jasmine tests verbosely, outputs info about individual tests
  grunt jasmine --filter=some          - runs jasmine tests for test specs matching 'some'
  grunt test                           - runs all tests
  grunt -v jasmine_nodejs:client_tests - run all client tests
  grunt -v jasmine_nodejs:server_tests - run all server tests

