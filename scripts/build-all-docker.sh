#!/usr/bin/bash

for dockerfilepath in $(find . -type f -name Dockerfile*);
do
    dockerfilename=$(basename $dockerfilepath) 
    pushd $(dirname $dockerfilepath)
    docker build -f $dockerfilename .
    if [ $? -ne 0 ];
    then
        echo -e "\e[1mBuilding $dockerfilepath failed\e[0m";
        exit 1
    fi 
    popd
done
