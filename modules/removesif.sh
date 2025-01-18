#!/bin/bash

for i in $(ls)
do
    if [ $i = "common" ] || [ $i == "test" ] || [ $i == "tmp" ] ; then
	continue
    fi
    if [ -d $i ] ; then
        if [ -f $i/$i.sif ] ; then
            rm -f $i/$i.sif
        fi
    fi
done
