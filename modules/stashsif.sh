#!/bin/bash

for i in $(ls)
do
    if [ $i == "test" ] || [ $i == "tmp" ] ; then
        continue
    fi
    if [ -d $i ] ; then
        for j in $(ls $i/*.sif 2> /dev/null)
        do
            echo "stashing $j"
            mv -f $j $j.stash
        done
    fi
done
