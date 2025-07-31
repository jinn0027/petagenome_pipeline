#!/bin/bash

NEED_MODULES=$(pwd)/../etc/need_modules.txt

function check_module () {
    module=$1
    if [ ! -f ${NEED_MODULES} ] || grep -w -q ${module} ${NEED_MODULES} ; then
        return 0
    else
        return 1
    fi
}

rm -f makeall.log

if groups | grep -q wheel ; then
    # SElinuxをpermissiveに設定。
    # これをやらないとapptainerでfakerootでdnfをやろうとするとエラーが発生することがある。
    if [ "$(getenforce)" != "Permissive" ] && [ "$(getenforce)" != "Disabled" ] ; then
        sudo setenforce 0
    fi
fi

make -j$(nproc) -C common all

for i in $(ls)
do
    if [ $i = "common" ] || [ $i == "test" ] || [ $i == "tmp" ] ; then
        continue
    fi
    
    if ! check_module $i ; then
	continue
    fi
    
    if [ -d $i ] ; then
        echo "$(date) $i" | tee -a makeall.log
	    make -j$(nproc) -C $i all
    fi
done
echo "$(date) complete" | tee -a makeall.log
