#!/bin/bash

rm -f makeall.log

# SElinuxをpermissiveに設定。
# これをやらないとapptainerでfakerootでdnfをやろうとするとエラーが発生する。
if [ "$(getenforce)" != "Permissive" ] && [ "$(getenforce)" != "Disabled" ] ; then
    sudo setenforce 0
fi

make -j$(nproc) -C common all

for i in $(ls)
do
    if [ $i = "common" ] || [ $i == "test" ] || [ $i == "tmp" ] ; then
        continue
    fi
    if [ -d $i ] ; then
        echo "$(date) $i" | tee -a makeall.log
	    make -j$(nproc) -C $i all
    fi
done
echo "$(date) complete" | tee -a makeall.log
