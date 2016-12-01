#!/usr/bin/env perl

use strict;
use warnings;

local @ARGV = ('httpd.conf');
local $^I = '.bak';
while (<>) {
    s{CUSTOM_PORT[\s]+80$}{CUSTOM_PORT  9090};
    s{^LoadModule auth_openidc}{\#LoadModule auth_openidc};
    s{^\#User(.*)$}{User www-data};
    s{^\#Group(.*)$}{Group www-data};
    s{^OIDC}{\#OIDC};

    s{^\#ServerName www.example.com:80}{ServerName localhost};

    s{^LoadModule ssl_module}{\#LoadModule ssl_module};
    s{^LoadModule proxy_connect_module}{\#LoadModule proxy_connect_module};
    s{^SSL}{\#SSL};

    s{^Define RSOCKET_PATH(.*)$}{Define NPG_DOCKER 1};
    next if (m{^<(.*?)Location(.*)}
          or m{^( *)<(.*?)Limit(.*)}
          or m{^( *)AuthType openid-connect}
          or m{^( *)Require valid-user});
    print;
}
