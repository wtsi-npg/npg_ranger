FROM ubuntu:14.04
MAINTAINER Andrew Nowak <an8@sanger.ac.uk>

# Hold versions in environment variables.
ENV NODE_VERSION 4.5.0
ENV MONGODB_VERSION 3.2.4
ENV HTSLIB_VERSION 1.3.1-npg-Apr2016
ENV SAMTOOLS_VERSION 1.3.1-npg-May2016
ENV BIOBAMBAM_VERSION 2.0.50-release-20160705161609
ENV FREEBAYES_VERSION 1.0.2-npg-Aug2016
ENV NPG_RANGER_VERSION 0.4.0

RUN apt-get update -qq \
  && apt-get install -qq \
    wget \
    build-essential \
    git \
    autoconf \
    zlib1g-dev \
    pkg-config \
    cmake

# Get Node.js
RUN wget "https://nodejs.org/dist/v$NODE_VERSION/node-v$NODE_VERSION-linux-x64.tar.xz" \
  && tar xJf "node-v$NODE_VERSION-linux-x64.tar.xz" -C /usr/local --strip-components=1 \
  && rm "node-v$NODE_VERSION-linux-x64.tar.xz"

# Get Mongodb
RUN wget "http://fastdl.mongodb.org/linux/mongodb-linux-x86_64-${MONGODB_VERSION}.tgz" \
  && tar xzf "mongodb-linux-x86_64-${MONGODB_VERSION}.tgz" -C /usr/local --strip-components=1 \
  && rm "mongodb-linux-x86_64-${MONGODB_VERSION}.tgz"

# Get HTSlib
RUN git clone --branch $HTSLIB_VERSION --depth 1 "https://github.com/wtsi-npg/htslib.git" htslib
WORKDIR ./htslib
RUN autoreconf -fi && ./configure --enable-plugins && make
WORKDIR ..

# Get SAMtools
RUN git clone --branch $SAMTOOLS_VERSION --depth 1 "https://github.com/wtsi-npg/samtools.git" samtools
WORKDIR ./samtools
RUN mkdir -p acinclude.m4 \
  && wget -o /dev/null -O acinclude.m4/ax_with_htslib.m4 "https://github.com/samtools/samtools/files/62424/ax_with_htslib.m4.txt" \
  && wget -o /dev/null -O acinclude.m4/ax_with_curses.m4 "http://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=m4/ax_with_curses.m4;hb=0351b066631215b4fdc3c672a8ef90b233687655" \
  && aclocal -I acinclude.m4 && autoreconf -fi \
  && ./configure --enable-plugins --without-curses \
  && make && make install
WORKDIR ..
RUN rm -rf samtools

# Get Biobambam
RUN wget "https://github.com/gt1/biobambam2/releases/download/${BIOBAMBAM_VERSION}/biobambam2-${BIOBAMBAM_VERSION}-x86_64-etch-linux-gnu.tar.gz" -O biobambam2.tar.gz \
  && tar xzf biobambam2.tar.gz -C /usr/local --strip-components=1 \
  && rm -r biobambam2.tar.gz

# Get Freebayes
RUN git clone --branch v$FREEBAYES_VERSION --depth 1 "https://github.com/wtsi-npg/freebayes.git" freebayes
WORKDIR freebayes
RUN git submodule update --init --recursive && make BAMTOOLS_ARGS="-DEnableNodeJS=true" && make install
WORKDIR ..
RUN rm -rf freebayes


# Get Ranger
RUN git clone --branch $NPG_RANGER_VERSION https://github.com/wtsi-npg/npg_ranger.git npg_ranger
WORKDIR npg_ranger
RUN npm install -g grunt-cli && npm install

