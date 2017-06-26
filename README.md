# npg_ranger project

[![Build Status](https://travis-ci.org/wtsi-npg/npg_ranger.svg?branch=master)](https://travis-ci.org/wtsi-npg/npg_ranger)
[![NPM dependencies](https://david-dm.org/wtsi-npg/npg_ranger.svg)](https://david-dm.org/wtsi-npg/npg_ranger)

## Introduction

This project is a contribution towards development of the streaming API under the
[GA4GH](http://ga4gh.org) initiative (GA4GH Directory API and Streaming). The project has
a server and a compatible command-line http client; both conform to [GA4GH API](https://github.com/samtools/hts-specs/blob/master/htsget.md)

## Package structure

```
 bin - node.js scripts, server and client
 lib - application libraries
 test - tests
 docker - scripts for deploying a server in docker containers in a cloud 
 docs - technical manuals
 docs/api - generated API (run "grunt doc" to generate)
 package.json - npm package definition
 Gruntfile.js - grunt runner specification
```

# Related projects

 - OAuth2 authorization server https://github.com/wtsi-npg/npg_sentry
 - GA4GH-compatible genome browser (biodalliance) extension https://github.com/wtsi-npg/dalliance