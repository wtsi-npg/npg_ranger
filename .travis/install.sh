#!/bin/bash

set -e -x

pushd /tmp
# mongo
wget "http://fastdl.mongodb.org/linux/mongodb-linux-x86_64-${MONGODB_VERSION}.tgz"
tar xfz "mongodb-linux-x86_64-${MONGODB_VERSION}.tgz"

# conda
mkdir -p software
wget -q "https://repo.continuum.io/miniconda/Miniconda2-${MINICONDA_VERSION}-Linux-x86_64.sh" -O "miniconda-${MINICONDA_VERSION}.sh"
/bin/sh "miniconda-${MINICONDA_VERSION}.sh" -b -p "${MINICONDA_HOME}"


# htslib & samtools
if [ ! "$(ls -A htslib)" ]; then
git clone --branch "${HTSLIB_VERSION}" --depth 1 https://github.com/samtools/htslib.git htslib
pushd htslib
autoreconf -fi
git reset --hard
./configure --prefix=/tmp/local --enable-libcurl
make
make install
popd
fi

if [ ! "$(ls -A samtools)" ]; then
git clone --branch "${SAMTOOLS1_VERSION}" --depth 1 https://github.com/samtools/samtools.git samtools
pushd samtools
mkdir -p acinclude.m4
pushd acinclude.m4
curl -L https://github.com/samtools/samtools/files/62424/ax_with_htslib.m4.txt > ax_with_htslib.m4
curl -L 'http://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=m4/ax_with_curses.m4;hb=0351b066631215b4fdc3c672a8ef90b233687655' > ax_with_curses.m4
popd
aclocal -I acinclude.m4
autoreconf -i
git reset --hard
LIBS='-lcurl -lcrypto -lssl' ./configure --prefix=/tmp/local --with-htslib=/tmp/htslib --without-curses
make
popd
fi

# biobambam
wget "https://github.com/gt1/biobambam2/releases/download/${BIOBAMBAM_VERSION}/biobambam2-${BIOBAMBAM_VERSION}-x86_64-etch-linux-gnu.tar.gz" -O biobambam2.tar.gz
mkdir biobambam2
tar xzf biobambam2.tar.gz -C biobambam2 --strip-components 1

# freebayes
if [ ! "$(ls -A freebayes)" ]; then
git clone --branch "${FREEBAYES_VERSION}" --depth 1 https://github.com/ekg/freebayes.git freebayes
pushd freebayes
git submodule update --init --recursive
make
popd
fi

# symlink to path
mkdir -p /tmp/usr/bin

ln -s "/tmp/mongodb-linux-x86_64-${MONGODB_VERSION}/bin/mongo" /tmp/usr/bin/mongo
ln -s "/tmp/mongodb-linux-x86_64-${MONGODB_VERSION}/bin/mongod" /tmp/usr/bin/mongod
ln -s "/tmp/mongodb-linux-x86_64-${MONGODB_VERSION}/bin/mongoimport" /tmp/usr/bin/mongoimport
ln -s "/tmp/miniconda/bin/conda" /tmp/usr/bin/conda # ? Surely since it was installed into software folder before this all might change too?
ln -s /tmp/samtools/samtools /tmp/usr/bin/samtools
ln -s /tmp/biobambam2/bin/bamstreamingmarkduplicates /tmp/usr/bin/bamstreamingmarkduplicates
ln -s /tmp/biobambam2/bin/bamseqchksum /tmp/usr/bin/bamseqchksum
ln -s /tmp/freebayes/bin/freebayes /tmp/usr/bin/freebayes

popd
