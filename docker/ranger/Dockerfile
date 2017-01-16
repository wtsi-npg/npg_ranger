FROM ubuntu:14.04
MAINTAINER Andrew Nowak <an8@sanger.ac.uk>

# Hold versions in environment variables.
ENV NODE_VERSION=4.5.0 \
    MONGODB_VERSION=3.2.4 \
    HTSLIB_VERSION=develop \
    SAMTOOLS_VERSION=develop \
    BIOBAMBAM_VERSION=2.0.50-release-20160705161609 \
    FREEBAYES_VERSION=1.1.0 \
    NPG_RANGER_VERSION=1.0.0 \
    APT_DEPENDENCIES="wget build-essential git autoconf pkg-config cmake"

RUN apt-get update -qq \
  && apt-get install -qq zlib1g-dev libcurl4-openssl-dev libssl-dev libbz2-dev liblzma-dev s3cmd curl \
  && apt-get update -qq \
  && apt-get install -qq $APT_DEPENDENCIES \

# Get Node.js
  && wget -q "https://nodejs.org/dist/v$NODE_VERSION/node-v$NODE_VERSION-linux-x64.tar.xz" \
  && tar xJf "node-v$NODE_VERSION-linux-x64.tar.xz" -C /usr/local --strip-components=1 \
  && rm "node-v$NODE_VERSION-linux-x64.tar.xz" \

# Get Mongodb - only needed for running grunt tests
# Use a dedicated mongodb container to host db in production
  && wget -q "http://fastdl.mongodb.org/linux/mongodb-linux-x86_64-${MONGODB_VERSION}.tgz" \
  && tar xzf "mongodb-linux-x86_64-${MONGODB_VERSION}.tgz" -C /usr/local --strip-components=1 \
  && rm "mongodb-linux-x86_64-${MONGODB_VERSION}.tgz" \

# Get HTSlib
  && git clone --branch $HTSLIB_VERSION "https://github.com/samtools/htslib.git" htslib \
  && cd ./htslib \
  # && git checkout c8475e824906a914f3aeeacef89f555009af9f2b \
  && autoreconf -fi \
  && git reset --hard c8475e824906a914f3aeeacef89f555009af9f2b \
  && ./configure --enable-libcurl \
  && make && make install \
  && cd .. \

# Get SAMtools
  && git clone --branch $SAMTOOLS_VERSION --depth 1 "https://github.com/samtools/samtools.git" samtools \
  && cd ./samtools \
  && mkdir -p acinclude.m4 \
  && wget -q -o /dev/null -O acinclude.m4/ax_with_htslib.m4 "https://github.com/samtools/samtools/files/62424/ax_with_htslib.m4.txt" \
  && wget -q -o /dev/null -O acinclude.m4/ax_with_curses.m4 "http://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=m4/ax_with_curses.m4;hb=0351b066631215b4fdc3c672a8ef90b233687655" \
  && aclocal -I acinclude.m4 && autoreconf -fi \
  && git reset --hard \
  && LIBS='-lcurl -lcrypto -lssl' ./configure --with-htslib=/htslib --without-curses \
  && make && make install \
  && cd .. \
  && rm -rf samtools \
  
# Get Biobambam
  && wget -q "https://github.com/gt1/biobambam2/releases/download/${BIOBAMBAM_VERSION}/biobambam2-${BIOBAMBAM_VERSION}-x86_64-etch-linux-gnu.tar.gz" -O biobambam2.tar.gz \
  && tar xzf biobambam2.tar.gz -C /usr/local --strip-components=1 \
  && rm -r biobambam2.tar.gz \

# Get Freebayes
  && git clone --branch v$FREEBAYES_VERSION --depth 1 "https://github.com/ekg/freebayes.git" freebayes \
  && cd ./freebayes \
  && git submodule update --init --recursive \
  && make \
  && make install \
  && cd .. \
  
# Get Ranger
  && npm install -g npg_ranger@${NPG_RANGER_VERSION}

COPY ./certs/*.pem /usr/share/ca-certificates/ranger/
RUN  if [ -d "/usr/share/ca-certificates/ranger" ]; \
  then \
    for i in /usr/share/ca-certificates/ranger/*.pem; \
    do \
      echo "$i"; \
      echo "ranger/$(basename $i)" | tee -a /etc/ca-certificates.conf; \
    done; \
    update-ca-certificates -v; \
  fi

WORKDIR /npg_ranger

RUN adduser --home /home/apprunner --disabled-password --gecos "" apprunner
COPY ./s3cfg /home/apprunner/.s3cfg
RUN chown -R apprunner:apprunner /npg_ranger /home/apprunner

# Remove build-time dependencies
RUN apt-get remove -qq $APT_DEPENDENCIES \
#  && apt-get autoremove --purge -qq \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

ADD https://raw.githubusercontent.com/vishnubob/wait-for-it/master/wait-for-it.sh /wait-for-it.sh

RUN chmod 0555 /wait-for-it.sh
EXPOSE 5050

USER apprunner

CMD [ "/bin/bash" ]

