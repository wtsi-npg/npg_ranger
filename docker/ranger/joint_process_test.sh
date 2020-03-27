#!/usr/bin/env sh

set -e

/expect_script_test.sh || echo "failed expect_script_test"

npg_ranger_server "$@"
