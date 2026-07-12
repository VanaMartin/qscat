#!/bin/bash

if [ -d "output" ]; then
    echo "Output directory exists"
    sdir=""
    if [ "$#" == 1 ]; then
        if [ -d "output/$1" ]; then 
            echo "Subdirectory $1 exists";
            sdir="$1"
        else 
            echo "Subdirectory missing"
            exit 1
        fi 
    fi 

    find output/$sdir -type f
    echo "Delete these files (y/n):"
    read answ
    if [ $answ=="y" ]; then
        find output/$sdir -type f -delete
    fi
fi
