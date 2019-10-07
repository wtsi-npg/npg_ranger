#!/bin/bash

set -e -x

pushd /tmp
# mongo
wget "http://fastdl.mongodb.org/linux/mongodb-linux-x86_64-${MONGODB_VERSION}.tgz"
tar xfz "mongodb-linux-x86_64-${MONGODB_VERSION}.tgz"

# conda
mkdir -p "${SOFTWARE_HOME}"
wget -q "https://repo.continuum.io/miniconda/Miniconda2-${MINICONDA_VERSION}-Linux-x86_64.sh" -O "miniconda-${MINICONDA_VERSION}.sh"
/bin/sh "miniconda-${MINICONDA_VERSION}.sh" -b -p "${MINICONDA_HOME}"
export PATH="${MINICONDA_HOME}/bin:$PATH"
export CONDA_ALWAYS_YES="true"

# samtools w/ conda
if [ ! "$(ls -A samtools)" ]; then
conda create -p "${SOFTWARE_HOME}/samtools/${SAMTOOLS1_VERSION}"
conda install -p "${SOFTWARE_HOME}/samtools/${SAMTOOLS1_VERSION}" -c "${CONDA_GENERIC_CHANNEL}" samtools="${SAMTOOLS1_VERSION}"
export PATH="${SOFTWARE_HOME}/samtools/${SAMTOOLS1_VERSION}/bin:$PATH"
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
ln -s "${SOFTWARE_HOME}"/samtools/"${SAMTOOLS1_VERSION}"/bin/samtools /tmp/usr/bin/samtools
#ln -s /tmp/samtools/samtools /tmp/usr/bin/samtools
ln -s /tmp/biobambam2/bin/bamstreamingmarkduplicates /tmp/usr/bin/bamstreamingmarkduplicates
ln -s /tmp/biobambam2/bin/bamseqchksum /tmp/usr/bin/bamseqchksum
ln -s /tmp/freebayes/bin/freebayes /tmp/usr/bin/freebayes

popd
