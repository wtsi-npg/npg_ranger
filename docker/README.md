# Docker-compose containers
This directory contains a collection of Dockerfiles, scripts, and a
docker-compose.yml which will allow you to quickly set up a working ranger
server.

This configuration is working with npg_ranger v2.0.0 in the ranger container. It
includes proof of concept Dockerfiles and docker-compose.yml to bring up a
mongodb container, an npg_ranger container and an apache container to act as a
reverse proxy.

## Preparing the instance
Instance preparation may be covered by an existing image. Otherwise these steps
should prepare the instance for the rest of the deployment.

```
export DEBIAN_FRONTEND=noninteractive
sudo apt-get update

sudo apt-get install htop s3cmd

# To install docker packages signed with the key
sudo apt-get -y install apt-transport-https ca-certificates
sudo apt-key adv --keyserver hkp://ha.pool.sks-keyservers.net:80 --recv-keys 58118E89F3A912897C070ADBF76221572C52609D

# cat /etc/issue to find ubuntu version for the correct ubuntu version in this case "ubuntu-xenial"
echo "deb https://apt.dockerproject.org/repo ubuntu-xenial main" | sudo tee /etc/apt/sources.list.d/docker.list
sudo apt-get update
sudo apt-get -y install linux-image-extra-$(uname -r) linux-image-extra-virtual
sudo apt-get -y install docker-engine
sudo apt-get -y upgrade

# add the user (ubuntu) to the docker group so sudo is not required for every docker command
sudo usermod -a -G docker $USER

sudo service docker start

# get docker compose binary
sudo curl -L https://github.com/docker/compose/releases/download/1.9.0/docker-compose-`uname -s`-`uname -m` -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose
```

## References volume

The deployment instructions assume there will be a **reference** root path to
pass to the containers. If an existing volume with references is available then
it may just need to be mounted.

```
# Find volume with references and attach to instance.
# Then mount in the instance:
sudo mkdir /references
sudo chmod a+r /references
# run lsblk to find references volume name then mount with something like:
sudo mount /dev/xvdf /references
```

You may want to save these changes to the fstab file so they don't need to be
repeated after stopping/starting the instance.

**Warning:** Be specially careful when executing this instruction as it will
prevent starting the instance in the future if wrong configuration is used.

```
# this line matches previous configuration, change if needed
sudo echo -e "/dev/xvdf\t/references\text4\tdefaults\t0 2" | sudo tee --append /etc/fstab
```
## Data volume

Similar to the case of the references, a data mounting point will be
required as part of the configuration for the containers. Here is where
db files will be saved.

```
# Find data volume attach to instance.
# Then mount in the instance:
sudo mkdir /data
# run lsblk to find references volume name then mount with something like:
sudo mount /dev/xvdg /data
```
Save changes to fstab if needed.

## Configuration for s3cmd
If data files are sourced from S3 (as in the sample data provided), s3cmd
configuration needs to be set.

For **AWS** use s3cmd interactive configruation

```
# interactive configuration to save s3cmd configuration at home
s3cmd --configure

```

For **OpenStack** use specific configuration file with internal entry points

## Get a clone of the project

```
# use specific branch if needed
git clone -b devel https://github.com/wtsi-npg/npg_ranger.git npg_ranger && pushd npg_ranger/docker
```

## Place configuration files in expected locations

### Configuration for s3cmd
Place the s3cmd configuration file generated in previous steps where it can
be used by docker. During container building phase this file will be passed to
the container.

```
cp ~/.s3cfg ranger/s3cfg
```

### Base configuration for apache
Place the httpd.conf in path where docker can see it.
```
cp ../docs/apache/httpd.conf rangerproxy/httpd.conf
```

### Working with temporary DNS hostname/IP
If the instance gets a new IP/DNS hostname regularly set environment variables for containers
to be aware if the current IP/hostname.

#### Using AWS meta data
Every time the intance is restarted and/or DNS/IP configuration changes

```
cp ranger/config_aws.json ranger/config.json # Get original configuration file

export PUBLIC_DNS_HOSTNAME=$(curl -s http://169.254.169.254/latest/meta-data/public-hostname) # To get public DNS hostname
export PUBLIC_IP=$(curl -s http://169.254.169.254/latest/meta-data/public-ipv4) # To get public IP

sed -i "s/<%PUBLIC_DNS_HOSTNAME%>/$PUBLIC_DNS_HOSTNAME/g" ranger/config.json # Set dns hostname in configuration file
sed -i "s/<%PUBLIC_IP%>/$PUBLIC_IP/g" ranger/config.json # Set public ip in configuration file
```

#### Openstack
```
cp ranger/config_os.json ranger/config.json # Get original configuration file

export PUBLIC_IP=<floating IP> # The floating IP of the machine
sed -i "s/<%PUBLIC_IP%>/$PUBLIC_IP/g" ranger/config.json # Set public ip in configuration file
```

### Add private CA certificates to ranger container
Place CA certificates in individual **.pem** files at *ranger/certs*. The
certificates will be copied into the ranger container and installed while it
is being build.
```
# e.g.
cp /usr/share/ca-certificates/local/my_cert.pem ranger/certs
```

## Use docker-compose
Use docker-compose to build and bring up containers.

In the generic form (AWS)
```
docker-compose up -d --build
```

Or supplying file for specific platform (OpenStack):
```
# If needed, edit docker-compose_os.yml to add IPs for DNS servers
docker-compose --file docker-compose_os.yml up -d --build
```

## Test data
A sample data file is provided. It can be used to test the service but needs to
be loaded.
```
# copy sample data to data mounting point
cp ./rangerdb/fileinfo.json /data/mongo/fileinfo.json

# Load some data (run only once! multiple times will create duplicate rows in database)
docker exec docker_rangerdb_1 mongoimport --db imetacache --collection fileinfo --jsonArray --file /data/db/fileinfo.json --upsert --upsertFields data_object,filepath_by_host
# This syntax will change in mongo 3.4 to something like --mode upsert --upsertFields data_object,filepath_by_host
```

## Testing deploy with sample data

```
echo "03555f613ce1c9cba69a862137f13b76  temp.bam" >> test_data.md5
echo "70b8f5f160e27210169ff50b013a75bb  temp.cram" >> test_data.md5
echo "f03053c6a581f18f1eb41a259e31a9b0  temp.sam" >> test_data.md5

curl "http://localhost:9090/npg_ranger/file?name=20818_1%23888.bam" -o temp.bam
curl "http://localhost:9090/npg_ranger/sample?accession=NA12878&format=cram&region=chr22:16100000-16105000" -o temp.cram
curl "http://localhost:9090/npg_ranger/sample?accession=NA12878&format=sam&region=chr22:16100000-16105000" -o temp.sam

md5sum -c test_data.md5 && rm temp.bam temp.cram temp.sam test_data.md5
```

These dockerfiles are a proof of concept only; no security has been enabled on
them, and are obviously not production-ready.
