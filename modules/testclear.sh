#!/bin/bash

for i in $(ls)
do
    if [ $i = "common" ] || [ $i == "test" ] || [ $i == "tmp" ] ; then
	continue
    fi
    if [ -f $i/test/c.sh ] ; then
	pushd $i/test >& /dev/null
	./c.sh
	popd >& /dev/null
    fi
done
