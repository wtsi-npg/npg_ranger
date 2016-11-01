#!/usr/bin/env perl

use strict;
use warnings;

local @ARGV = ('httpd.conf');
local $^I = '.bak';
while (<>) {
    s{CUSTOM_PORT  80$}{CUSTOM_PORT  9090};
    s{^LoadModule auth_openidc}{\#LoadModule auth_openidc};
    s{^\#User(.*)$}{User www-data};
    s{^\#Group(.*)$}{Group www-data};
    s{^OIDC}{\#OIDC};
    s{^Define RSOCKET_PATH(.*)$}{Define NPG_DOCKER 1};
    next if (m{^<(.*?)Location(.*)}
          or m{^( *)<(.*?)Limit(.*)}
          or m{^( *)AuthType openid-connect}
          or m{^( *)Require valid-user});
    print;
}
