#!/bin/bash

rm -f makesbx.log
for i in $(ls)
do
    if [ $i = "common" ] || [ $i == "test" ] || [ $i == "tmp" ] ; then
	continue
    fi
    if [ -d $i ] ; then
	    #make -C $i all
        echo "$(date) $i" | tee -a makesbx.log
        make -C $i ${i}.sbx
    fi
done
echo "$(date) complete" | tee -a makesbx.log
