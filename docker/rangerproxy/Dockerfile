FROM httpd:2.4.23
MAINTAINER Andrew Nowak <an8@sanger.ac.uk>

ENV NPG_RANGER_VERSION=1.0.0

ADD https://raw.githubusercontent.com/jmtcsngr/npg_ranger/docker_update/docs/apache/httpd.conf /usr/local/apache2/httpd.conf
COPY ./docker_conf.pl /usr/local/apache2/docker_conf.pl

RUN perl docker_conf.pl \
  && mv /usr/local/apache2/httpd.conf /usr/local/apache2/conf/httpd.conf

EXPOSE 9090
