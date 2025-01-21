#!/bin/bash

fq1=../../test/s_6_1.fastq.gz
fq2=../../test/s_6_2.fastq.gz

odir=results

refdir=ref

log=t.log

ret=0

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} ../bbmap.sif bbmap -h > ${log} 2>&1
for i in $(ls $odir/*.html)
do
    j=$refdir/$(basename $i)
    diff -q $i $j >> ${log} 2>&1 && :
    if [ $? -ne 0 ]; then
        ret=1
    fi
done

if [ ${ret} -eq 0 ]; then
    echo " PASSED" | tee -a ${log}
else
    echo " FAILED" | tee -a ${log}
fi

exit ${ret}

# #{bbmap_dir}/bbmap/bbmap.sh -Xmx#{mem}g threads=#{n_threads} ref=#{contig_}
# #{bbmap_dir}/bbmap/bbmap.sh -Xmx#{mem}g threads=#{n_threads} in=#{query1_} in2=#{query2_} ambiguous=random minid=0.95 pairlen=1500 out=#{out}.sam scafstats=#{out}.scafstats
#10_contig_coverage.rb













