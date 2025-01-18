#!/bin/bash

s=(
    virsorter
)

f=(
    virsorter.tar.gz
)

k=(
    virsorter
)

n=$(( ${#s[@]} - 1 ))

for i in $(seq 0 $n)
do
    #echo ${s[$i]}
    cp -r tmp ${s[$i]}
    pushd ${s[$i]}
    sed -i "s#SSS#${s[$i]}#g" *
    sed -i "s#FFF#${f[$i]}#g" *
    sed -i "s#KKK#${k[$i]}#g" *
    mv SSS.def ${s[$i]}.def
    mv SSS.sh ${s[$i]}.sh
    popd
done
