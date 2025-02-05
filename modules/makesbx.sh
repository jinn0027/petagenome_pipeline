#!/bin/bash

rm -f makesbx.log
for i in $(ls)
do
    if [ $i = "common" ] || [ $i == "test" ] || [ $i == "tmp" ] ; then
        continue
    fi
    if [ -d $i ] ; then
        echo "$(date) $i" | tee -a makesbx.log
        make -j$(nproc) -C $i ${i}.sbx
    fi
done
echo "$(date) complete" | tee -a makesbx.log
