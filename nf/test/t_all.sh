#!/bin/bash

for l in $(seq 1 3)
do
    for i in $(ls ../lv${l}/*.nf | sort)
    do
        j=$(basename $i)
        ./t.sh $j
    done
done

##############

tot=0
num_pass=0

F1=.l
F2=.m
for dir in $(ls ref)
do
    if [ ! -d ref/$dir ] ; then
        continue
    fi
    fail=0
    for i in $(find ref/$dir -type f -print | sort)
    do
        j=$(echo $i | sed "s#^ref/#out/#")
        echo ${i} | grep -e ".gz\$" >& /dev/null && :
        if [ $? -eq 0 ] ; then
            gunzip -c ${i} > ${F1}
            gunzip -c ${j} > ${F2} 
            diff -q ${F1} ${F2} >& /dev/null
        else
            echo ${i} | grep -e ".zip\$" >& /dev/null && :
            if [ $? -eq 0 ] ; then
                unzip -p ${i} > ${F1}
                unzip -p ${j} > ${F2}
                diff -q ${F1} ${F2} >& /dev/null
            else
                diff -q $i $j >& /dev/null
            fi
        fi
        if [ $? -ne 0 ] ; then
            echo $i
            exit 1
            fail=$(( fail + 1 ))
        fi
    done
    if [ ${fail} -eq 0 ] ; then
        echo "$dir PASS"
        num_pass=$(( num_pass + 1 ))
    else
        echo "$dir FAIL"
    fi
    tot=$(( tot + 1 ))
done

rm -f ${F1} ${F2}

echo "PASS=${num_pass} / ${tot}"


