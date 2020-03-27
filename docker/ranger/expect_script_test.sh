#!/usr/bin/expect -f
set timeout 10

spawn iinit

expect {
    "Enter your current iRODS password:" {
        send -- "$env(IRODS_PASSWORD)\r"
        exp_continue
    }
    "failed with error" {
        puts "Error in entering iRODS password"
        exit 1
        interact
    }
}
