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

