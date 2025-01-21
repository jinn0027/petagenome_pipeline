#!/bin/bash

n_threads=$(nproc)
mem=128

ref=../../test/ecoli_1K_1.fa.gz
fq1=../../test/s_6_1.fastq.gz
fq2=../../test/s_6_2.fastq.gz

wdir=work
odir=results

refdir=ref

log=t.log

ret=0

ref=$(cd $(dirname ${ref}) && pwd)/$(basename ${ref})
fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})

wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir} ${wdir}
/usr/local/bin/apptainer exec --bind ${ref},${fq1},${fq2},${odir},${wdir} ../bbmap.sif bbmap.sh -Xmx${mem}g threads=${n_threads} ref=${ref} path=${wdir} > ${log} 2>&1

/usr/local/bin/apptainer exec --bind ${ref},${fq1},${fq2},${odir},${wdir} ../bbmap.sif bbmap.sh -Xmx${mem}g threads=${n_threads} path=${wdir} in=${fq1} in2=${fq2} ambiguous=random minid=0.95 pairlen=1500 out=${odir}/out.sam scafstats=${odir}/out.scafstats >> ${log} 2>&1

for i in $(ls $odir/*.sam ${odir}/*.scafstats)
do
    j=${refdir}/$(basename $i)
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













