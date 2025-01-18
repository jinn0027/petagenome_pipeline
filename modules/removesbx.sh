#!/bin/bash

for i in $(ls)
do
    if [ $i = "common" ] || [ $i == "test" ] || [ $i == "tmp" ] ; then
	continue
    fi
    if [ -d $i ] ; then
        if [ -d $i/$i.sbx ] ; then
            sudo rm -rf $i/$i.sbx
        fi
    fi
done
