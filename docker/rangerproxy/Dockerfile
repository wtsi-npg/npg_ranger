FROM httpd:2.4.41
MAINTAINER Grisha Hushchyn <gh13@sanger.ac.uk>

ENV APT_DEPENDENCIES="wget git xz-utils autoconf pkg-config libcurl4-openssl-dev libssl-dev libjansson-dev build-essential libcjose-dev libpcre++-dev" \
    APT_DEPENDENCIES_KEEP="ca-certificates"

RUN apt-get update -qq \
  && apt-get install -qq $APT_DEPENDENCIES \
  && apt-get install -qq $APT_DEPENDENCIES_KEEP \
\
# Get mod-auth-openidc
  && git clone https://github.com/zmartzone/mod_auth_openidc.git --branch master auth_openidc \
\
# Build mod-auth-openidc
  && cd auth_openidc \
  && ./autogen.sh \
  && ./configure --with-apxs2=/usr/local/apache2/bin/apxs \
  && make \
  && make install \
  && cd ..

RUN apt-get remove -qq $APT_DEPENDENCIES \
#  && apt-get autoremove --purge -qq \
  && apt-get clean \
  && rm -rf ./auth_openidc/* \
  && rm -rf /var/lib/apt/lists/*

EXPOSE 9090
