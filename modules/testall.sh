#!/bin/bash

for i in $(ls)
do
    if [ $i = "common" ] || [ $i == "test" ] || [ $i == "tmp" ] ; then
	continue
    fi
    if [ -f $i/test/t.sh ] && [ -d $i/test/ref ] ; then
	pushd $i/test >& /dev/null
	echo -n "$i : "
	./t.sh
	popd >& /dev/null
    fi
done
