FROM httpd:2.4.27
MAINTAINER Andrew Nowak <an8@sanger.ac.uk>

ENV NODE_VERSION=6.10.3 \
    NPG_RANGER_VERSION=1.0.0 \
    APT_DEPENDENCIES="wget build-essential git autoconf pkg-config cmake"

RUN apt-get update -qq \
  && apt-get install -qq $APT_DEPENDENCIES \

# Get Node.js
  && wget -q "https://nodejs.org/dist/v$NODE_VERSION/node-v$NODE_VERSION-linux-x64.tar.xz" \
  && tar xJf "node-v$NODE_VERSION-linux-x64.tar.xz" -C /usr/local --strip-components=1 \
  && rm "node-v$NODE_VERSION-linux-x64.tar.xz" \
  && npm update -g npm

RUN git clone --branch npg_ranger_master --depth 1 https://github.com/wtsi-npg/dalliance.git dalliance \
  && cd dalliance \
  && npm install \
  && ./node_modules/.bin/gulp \
  && ln -s node_modules/npg_ranger/lib/uiclient/index.html . \
  && cd ../htdocs \
  && ln -s ../dalliance .

COPY ./docker_conf.pl /usr/local/apache2/docker_conf.pl
COPY ./httpd.conf /usr/local/apache2/httpd.conf

RUN perl docker_conf.pl \
  && mv /usr/local/apache2/httpd.conf /usr/local/apache2/conf/httpd.conf

RUN apt-get remove -qq $APT_DEPENDENCIES \
#  && apt-get autoremove --purge -qq \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

EXPOSE 9090
