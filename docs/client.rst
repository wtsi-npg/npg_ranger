#############
Client Manual
#############

Running the client
==================

When built from source
----------------------

With output directed to stdout (by default)

::

  $ bin/client.js "http://192.168.0.1:5050/resources/A01?referenceName=1&start=167856&end=173507&format=BAM"

With output writen to file

::

  $ bin/client.js "http://192.168.0.1:5050/resources/A01?referenceName=1&start=167856&end=173507&format=BAM" output.bam


When installed with npm
-----------------------

If installed locally


With output directed to stdout (by default)

::

  $ ./node_modules/npg_ranger/bin/client.js "http://192.168.0.1:5050/resources/A01?referenceName=1&start=167856&end=173507&format=BAM"

Or

::

  $ ./node_modules/.bin/npg_ranger_client "http://192.168.0.1:5050/resources/A01?referenceName=1&start=167856&end=173507&format=BAM"

With output writen to file

::

  $ ./node_modules/npg_ranger/bin/client.js "http://192.168.0.1:5050/resources/A01?referenceName=1&start=167856&end=173507&format=BAM" output.bam

Or

::

  $ ./node_modules/.bin/npg_ranger_client "http://192.168.0.1:5050/resources/A01?referenceName=1&start=167856&end=173507&format=BAM" output.bam

If installed globally
---------------------

With output directed to stdout (by default)

::

  $ npg_ranger_client "http://192.168.0.1:5050/resources/A01?referenceName=1&start=167856&end=173507&format=BAM"

With output writen to file

::

  $ npg_ranger_client "http://192.168.0.1:5050/resources/A01?referenceName=1&start=167856&end=173507&format=BAM" output.bam


Options
=======

--accept-trailers
-----------------

Request http trailers from the server, if implemented. This will help indicate if errors occurred during transfer, after the initial http status has been transmitted.

--loglevel <level>
------------------

Output all logging messages up to and including the given priority level.
Known levels:   

- debug
- info
- warn
- error [default]

--token_config <path>
---------------------

Provide the path to an npg_sentry authorisation token, inside a JSON file of the form:

::

  { "npg_sentry_token_bearer": "abcdefghijkl...xyz" }
  
--with_ca <path>
----------------

Path to a PEM formatted CA certificate to verify the npg_sentry server's identity.
