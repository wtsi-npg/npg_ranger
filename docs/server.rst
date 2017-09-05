#############
Server Manual
#############

General
=======

In a production setup the main entities on the server side are a node.js ranger server,
an Apache httpd server acting as a reverse proxy, behind which the ranger server runs,
and a MongoDB database that keeps metadata about location of sequencing data. All three
are capable of handling many thousands queries simultaneously. If the number of simultaneous
requests is not capped, the load on the underlying sequencing data storage might
become too high. A configuration for the httpd server, supplied in this project,
limits the number of simultaneously run queries. It can be adjusted either way
depending on the type of underlying storage.

Running
=======

1. Connect to iRODs if not already connected

2. Ensure you have bioinformatics tools in path

   2.1 samtools v1.5 or higher

   2.2 biobambam v2.0.50

   2.3 freebayes `v1.1.0
   <https://github.com/ekg/freebayes/tree/v1.1.0>`_

3. Create `configuration file
   <https://github.com/wtsi-npg/npg_ranger/blob/master/docs/config.json>`_
   with parameters needed e.g. mongo database url, path for reference root and
   port

4. Run server

Parameters
----------

Parameters can be used to configure the server. You can use -h to list
supported parameters.

::

  bin/server.js -h

Parameters and extended configuration can be defined in a configuration
file. An example configuration file can be found at
`docs/config.json <https://github.com/wtsi-npg/npg_ranger/blob/master/docs/config.json>`_.

::

  bin/server.js -c yourConfig.json


Providing essential configuration
---------------------------------

An essential parameter to start the server is the mongo database url. You can
set this parameter by creating a configuration file and passing it to the
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

If a reference root is required to resolve full reference paths from the entries
in the database, the root path should be provided in the startup configuration.

::

 # in a configuration file

 {
   "mongourl":   "mongodb://<url of mongo server>:<port>/imetacache",
   "references": "/path_to_ref_root/"
 }

 # as parameter
 bin/server.js -r "/path_to_ref_root/"

It is possible to use ssl connections to the mongo server (if the server is
configured to support them). Modify the mongodb url passing options for ssl
and certificate checking.

::
 
 If working with a self signed certificate in the mongo server

 "mongourl": "mongodb://<server>:<port>/imetacache?ssl=true&sslValidate=false"

 If using a certificate signed by a know CA

 "mongourl": "mongodb://<server>:<port>/imetacache?ssl=true"

If using a *CA* which is not installed in the system trust store, it is possible
to let node know that extra certificates can be found in a specific ``pem`` 
file by setting an enviroment variable, see `node cli documentation 
<https://nodejs.org/api/cli.html#cli_node_extra_ca_certs_file>`_.

::

 export NODE_EXTRA_CA_CERTS=/<path-to>/CA.pem

Running a secure server using HTTPS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run an https server at least a certificate and a private key must be provided
in PEM format. If the server private key was generated with a passphrase, the
passphrase must be provided as part of the configuration. Paths to the pem
files with the private key and the certificate can be passed as start up
options using ``--secure_key`` and ``--secure_cert`` or by configuration file
using ``secure_key`` and ``secure_cert`` entries. If a passphrase is needed, it
can only be provided in the configuration file under the ``secure_passphrase``
entry. For security reasons both .pem files and the server configuration file
must have proper access permissions, e.g. ``chmod 400 config.json``.

Secure connection to authorisation service
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If authorisation is being provided by an external service (**npg_sentry**),
extra configuration should be provided to establish a secure connection to the
service. At least an entry point *URL* must be provided with the ``--authurl``
option. If the listening server is expecting for **npg_ranger** to establish the
connection using pre-arranged certificates, ``--auth_cert`` and ``--auth_key``
options should be used. If the key for client side authentication is passphrase
protected, ``auth_key_passphrase`` should be set in the configuration file. The
passphrase for the client side key can only be set in the configuration file,
there is no matching command line option for it. Finally if the certificate of
authorisation service is signed by a private *CA*, **npg_ranger** will require
access to the *CA*'s certificate to validate the certificate presented by the
authorisation service. A path to the *CA* can be configured using the option
``--auth_ca``.


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
 npg_ranger_server

 #on a custom socket
 npg_ranger_server /tmp/my.sock

 #on a custom port 9447
 npg_ranger_server 9447

 #on a custom port and skip authentication
 npg_ranger_server -p PORT -s

 #changing time to wait before killing child processes
 npg_ranger_server -g SECONDS


EXAMPLES AND COMPATIBLE CLIENTS
===============================

available urls
==============

There are three different url paths recognised by the server:

::

 /file?name=$NAME[&directory=$DIR]
 /sample?accession=$ACCESSION[&format={BAM,SAM,CRAM,VCF}][&region=$REG]
 /ga4gh/sample/$ACCESSION[&format={BAM,SAM,CRAM,VCF}][&referenceName=$CHR&start=$STARTPOS&end=$ENDPOS]
 # $REG is in format <referenceName>:<startLoc>-<endLoc>

Each will provide a response in a different way:

/file will search the database for a file with matching name, then will stream that file to you. This url can only return one file, so if there is more than one file with $NAME, you will be prompted to also provide a directory $DIR. This url supports byte-range serving using the Content-Range header.

/sample will search for content files with given accession, merge them, then stream the file (or specified region) in BAM format (unless overridden).

/ga4gh/sample will provide a json response, mapping the url to a /sample url with the same accession and queries. The npg_ranger client and `our biodalliance fork`__ will automatically follow this redirect, curl and other http clients will not.

.. _Biodall: https://github.com/wtsi-npg/dalliance

__ Biodall_

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

The default output format is BAM. Use 'format' option with value either 'SAM' or
'BAM' or 'CRAM' to change the output format.

nodejs client (this project)
----------------------------
bin/client.js
A simple trailer header aware client that works with a socket server.

Biodalliance
------------
A custom npg_ranger track is added to the Biodalliance genome browser
https://github.com/wtsi-npg/dalliance/tree/npg_ranger_master

Authentication and authorisation
================================

Authentication should be done by a front server. It is expected that the
incoming request has X-Remote-User header set. The data will be served if the
remote user has 'read' permission for all files that have to be merged/served.

APACHE REVERSE PROXY
====================

Setting up the server
---------------------

::

 wget http://mirrors.ukfast.co.uk/sites/ftp.apache.org//httpd/httpd-2.4.27.tar.gz
 tar -xzvf httpd-2.4.27.tar.gz
 cd httpd-2.4.27
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

If a reverse proxy is set as an entry point for the application, the server will
need to be aware of the reverse proxy addresses and paths mapped. The list of
addresses and paths can be provided in the configuration file.

::

  {
    "proxylist": {
      "http://server:port": "http://server:port/mapped_path"
    }
  }

Example configuration entries for an Apache reverse proxy can be found bellow:

::

  ProxyPreserveHost On
  # to a local server listening on a unix socket, requires Apache v 2.4.7 at least
  ProxyPass /        unix:/path_to/my.socket|http://localhost/
  ProxyPassReverse / unix:/path_to/my.socket|http://localhost/
  # To a local URL
  #ProxyPass /        http://localhost:9030/
  #ProxyPassReverse / http://localhost:9030/

CORS headers
------------

If the server needs to provide data for browser clients, CORS headers may need
to be configured. A list of allowed origins can be passed as part of the
configuration file.

::

 {
   "originlist": ["http://one_origin.com", "http://other_origin.com"]
 }

If it is not possible to enumerate the origins to be allowed, the least secure
option of allowing all origins can be configured at server startup with the
--anyorigin option.


Filtering results
=================
Use query parameters to modify your results


/sample route
-------------

Required:
~~~~~~~~~
accession: Show files from given sample accession number::

  accession=ABC123456

Optional:
~~~~~~~~~

format: Return files, transforming data to given format::

  format={sam,bam,cram,vcf}

region: Return files, showing only the given genomic region::

  region=chr13:1700000-1700200

Filtering:
~~~~~~~~~~

WARNING:

The default filter values should work for the majority of cases. Not using the default values will dramatically increase the chance of errors occurring; either through attempts to access forbidden data, or by attempting to merge files with non-matching references.

Each filter can take value 'undef' to search for files where the attribute corresponding to given filter does not exist. Each filter can take an empty string as a value to search for files without querying by the attribute corresponding to that filter.

Each filter can be suffixed with '_not' to search for files with any value of the attribute *except* the given filter value. Giving value 'undef' to this form of the filter will return all files where the attribute corresponding to the filter exists, regardless of value.

Not specifying a filter in the query will filter by the default value if it exists, or otherwise will ignore that filter (the same as giving an empty string, above)

| For Example:
| ``target=1`` searches for files where target=1
| ``target=`` searches for all files, does not query by target
| ``target=undef`` searches for files where target has not been defined
| ``target_not=X`` returns the inverse of ``target=X``


+------------------+-----------+---------------+
| filter name      | default   | common values |
+==================+===========+===============+
| target           | 1         | 1, 0, library |
+------------------+-----------+---------------+
| manual_qc        | 1         | 1, 0          |
+------------------+-----------+---------------+
| alignment        | 1         | 1, 0          |
+------------------+-----------+---------------+
| alt_target       | n/a       | 1             |
+------------------+-----------+---------------+
| alt_process      | n/a       |               |
+------------------+-----------+---------------+
| alignment_filter | n/a       | phix,human,...|
+------------------+-----------+---------------+
