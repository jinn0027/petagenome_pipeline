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

for i in $(ls)
do
    if [ ! -d $i/$i.sbx ] ; then
        continue
    fi
    
    if [ $i = "common" ] || [ $i == "test" ] || [ $i == "tmp" ] ; then
        continue
    fi

    if ! check_module $i ; then
	continue
    fi
    
    if [ -f $i/test/c.sh ] ; then
        pushd $i/test >& /dev/null
        ./c.sh
        popd >& /dev/null
    fi
done
