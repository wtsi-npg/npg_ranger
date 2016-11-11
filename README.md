# npg_ranger project

This branch is locked at the 0.6.0 release. It includes proof of concept Dockerfiles and
docker-compose.yml to bring up a mongodb container, an npg_ranger container and an apache
container to act as a reverse proxy. Instructions to do this are included in
[docker/README.md](docker/README.md). Further work will be necessary to bring these files
up-to-date with new changes to npg_ranger and how we are configuring apache.

[![Build Status](https://travis-ci.org/wtsi-npg/npg_ranger.svg?branch=master)](https://travis-ci.org/wtsi-npg/npg_ranger)
[![NPM dependencies](https://david-dm.org/wtsi-npg/npg_ranger.svg)](https://david-dm.org/wtsi-npg/npg_ranger)

## Introduction

This project is a contribution towards development of the streaming API and servers under the
[GA4GH](http://ga4gh.org) initiative (GA4GH Directory API and Streaming). The project has
both a server and a compatible http client.

The server conforms to as yet unpublished at the time of writing GA4GH API. It uses
the http transfer encoding protocol to stream sequencing data.

## Package structure

```
 bin - node.js scripts, server and client
 lib - application libraries
 test - tests
 docs - technical manuals
 docs/api - generated API (run "grunt doc" to generate)
 package.json - npm package definition
 Gruntfile.js - grunt runner specification
```
