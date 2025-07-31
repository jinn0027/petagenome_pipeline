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

rm -f testall.log
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
    
    echo "$(date) $i" | tee -a testall.log
    if [ -f $i/test/skip_test ] ; then
        echo "$i : SKIPPED" | tee -a testall.log
        continue
    fi
    if [ -f $i/test/t.sh ] && [ -d $i/test/ref ] ; then
        pushd $i/test >& /dev/null
        echo -n "$i : " | tee -a ../../testall.log
        ./t.sh | tee -a ../../testall.log
        popd >& /dev/null
    fi
done
echo "$(date) complete" | tee -a testall.log
