#!/bin/bash

for i in $(ls)
do
    if [ $i == "test" ] || [ $i == "tmp" ] ; then
        continue
    fi
    if [ -d $i ] ; then
        rm -f *.stash
    fi
done
