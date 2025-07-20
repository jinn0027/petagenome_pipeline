#!/bin/bash

for i in $(ls)
do
    if [ $i == "test" ] || [ $i == "tmp" ] ; then
        continue
    fi
    if [ -d $i ] ; then
        for j in $(ls $i/*.sif.stash 2> /dev/null)
        do
            k=${j%%.stash}
            echo "unstashing $k"
            mv -f $j $k
        done
    fi
done
