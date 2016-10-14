#############
Server Manual
#############

Running
=======

1. Connect to iRODs if not already connected

2. Ensure you have bioinformatics tools in path

   2.1 samtools v1.3 or higher

   2.2 biobambam v2.0.50

   2.3 freebayes `v1.0.2-npg-Aug2016 <https://github.com/wtsi-npg/freebayes/tree/v1.0.2-npg-Aug2016>`_

3. Create configuration file with parameters needed e.g. mongo database url, path for reference root and port

4. Run server

Parameters
----------

Parameters can be used to configure the server. You can use -h to list
supported parameters.

::

  bin/server.js -h

Parameters and extended configuration can be defined in a configuration
file. An example configuration file can be found at docs/config.json.

::

  bin/server.js -c yourConfig.json


Providing essential configuration
---------------------------------

An essential parameter to start the server is the mongo database url. You can set
this parameter by creating a configuration file and passing it to the
server.

::

 {
   "mongourl": "mongodb://<url of mongo server>:<port>/imetacache"
 }

 # from  a configuration file
 bin/server.js -c <...>/yourConfig.json

Or by passing the parameter when starting the server with the -m option.

::

 #providing url of mongo server
 bin/server.js -m 'mongodb://<url>:<port>/imetacache'

 #providing path to unix socket
 bin/server.js -m 'mongodb:///tmp/mongodb-27017.sock/imetacache'


Other options
-------------

Installed from source

::


 #on the default socket /tmp/${USER}/npg_ranger.sock
 bin/server.js

 #on a custom socket
 bin/server.js /tmp/my.sock

 #on a custom port 9447
 bin/server.js 9447

 #on a custom port and skip authentication
 bin/server.js -p PORT -s

 #changing time to wait before killing child processes
 bin/server.js -g SECONDS

Installed with npm

If installed locally

::

 #on the default socket /tmp/${USER}/npg_ranger.sock
 ./node_modules/.bin/npg_ranger_server

 #on a custom socket
 ./node_modules/.bin/npg_ranger_server /tmp/my.sock

 #on a custom port 9447
 ./node_modules/.bin/npg_ranger_server 9447

 #on a custom port and skip authentication
 ./node_modules/.bin/npg_ranger_server -p PORT -s

 #changing time to wait before killing child processes
 ./node_modules/.bin/npg_ranger_server.js -g SECONDS

If installed globally

::

 #on the default socket /tmp/${USER}/npg_ranger.sock
 npg_ranger_server.js

 #on a custom socket
 npg_ranger_server.js /tmp/my.sock

 #on a custom port 9447
 npg_ranger_server.js 9447

 #on a custom port and skip authentication
 npg_ranger_server.js -p PORT -s

 #changing time to wait before killing child processes
 npg_ranger_server.js -g SECONDS


EXAMPLES AND COMPATIBLE CLIENTS
===============================

curl
====

::

 curl -H "Content-type: application/octet-stream" -X "GET" 'localhost:9444/sample?region=Zv9_scaffold3541&accession=ERS1023809'

no files found - an empty reply
one file found - an outcome of samtools view
multiple files found - an outcome of samtools merge

::

 curl -H "Content-type: application/octet-stream" -X "GET" 'localhost:9444/file?directory=/seq/18691&region=Zv9_scaffold3541&irods=1&name=18691_1%231.cram'
 curl -H "Content-type: application/octet-stream" -X "GET" 'localhost:9444/file?directory=/staging/path&region=Zv9_scaffold3541&name=18691_1%231.cram'

The default output format is bam. Use 'format' option with value either 'SAM' or 'BAM' or 'CRAM' to change the output format.

nodejs client (this project)
----------------------------
bin/client.js
A simple trailer header aware client that works with a socket server.

Biodalliance
------------
A custom npg_ranger track is added to the Biodalliance genome browser
https://github.com/wtsi-npg/dalliance

Authentication and authorisation
================================

Authentication should be done by a front server. It is expected that the incoming request has X-Remote-User header set. The data will be served if the remote user has 'read' permission for alll files that have to be merged/served.

APACHE REVERSE PROXY
====================

Setting up the server
---------------------

::

 wget http://mirrors.ukfast.co.uk/sites/ftp.apache.org//httpd/httpd-2.4.18.tar.gz
 tar -xzvf httpd-2.4.18.tar.gz
 cd httpd-2.4.18
 ./configure --enable-load-all-modules --prefix=${HOME}/apache_build
 make
 make install
 cd ${HOME}/apache_build
 vi conf/httpd.conf # edit the file

 #from anywhere
 ${HOME}/apache_build/bin/httpd -k start

LDAP authorisation config
-------------------------

::

 # Limit authentication to only GET and POST requests, auth will not be sent with OPTIONS
 <Location / >
   <Limit GET POST>
     AuthType Basic
     AuthBasicProvider ldap
     AuthName "LDAP Login For NPG Streaming"
     AuthLDAPURL "sanger ldap string"
     Require valid-user
     AuthLDAPRemoteUserAttribute uid
     RewriteEngine On
     RewriteRule .* - [E=PROXY_USER:%{LA-U:REMOTE_USER},NS]
     RequestHeader set X-Remote-User %{PROXY_USER}e
   </Limit>
  </Location>

Reverse proxy configuration
---------------------------

::

  ProxyPreserveHost On
  # to a local server listening on a unix socket, requires Apache v 2.4.7 at least
  ProxyPass /        unix:/path_to/my.socket|http://localhost/
  ProxyPassReverse / unix:/path_to/my.socket|http://localhost/
  # To a local URL
  #ProxyPass /        http://localhost:9030/
  #ProxyPassReverse / http://localhost:9030/

  # Additional headers to forward
  RewriteEngine On
  RewriteCond "%{HTTPS}" =off
  RewriteRule ^\/npg_ranger\/.* - [E=XPROTOCOL:http]
  RewriteCond "%{HTTPS}" =on
  RewriteRule ^\/npg_ranger\/.* - [E=XPROTOCOL:https]
  # Use ":" as suffix of protocol eg "http:"
  RequestHeader set X-Forwarded-Proto  "%{XPROTOCOL}e:"

CORS headers
------------

::

 Header set Access-Control-Allow-Origin "SOME_SERVER_URL"
 Header set Access-Control-Allow-Methods "GET"
 Header set Access-Control-Allow-Credentials "true"

Or, if no authentication is necessary,

::

 Header set Access-Control-Allow-Origin "*"
 Header set Access-Control-Allow-Methods "GET"
