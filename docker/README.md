# Docker-compose containers
This directory contains a collection of Dockerfiles, scripts, and a docker-compose.yml which will allow you to quickly set up a working server and database behind an apache reverse proxy.

Simply `cd` to this directory, and then run
`docker-compose up`.

The server should now be running on localhost:9090, although the database is empty.

To populate the database, using an example fileinfo.json which provides metadata for the test data provided with the main server as an example:

```
docker-compose up -d
cp ./rangerdb/fileinfo.json ./rangerdb/volume/
docker exec docker_rangerdb_1 mongoimport --db imetacache --collection fileinfo --jsonArray --file /data/db/fileinfo.json
```

Now either run the ranger client or use your browser to visit
```
http://localhost:9090/npg_ranger/file?name=20818_1%23888.bam&format=sam
```
and a sam file should be returned.

These dockerfiles are a proof of concept only; no security has been enabled on them, and are obviously not production-ready.
