#!/bin/bash

pass=0
fail=0
for i in $(ls)
do
    if [ $i = "common" ] || [ $i == "test" ] || [ $i == "tmp" ] ; then
	continue
    fi
    if [ -d $i ] ; then
        if [ -f $i/$i.sif ] ; then
            echo "o $i"
            pass=$(( pass + 1 ))
        else
            echo "x $i"
            fail=$(( fail + 1 ))
        fi
    fi
done

echo "pass : $pass"
echo "fail : $fail"
