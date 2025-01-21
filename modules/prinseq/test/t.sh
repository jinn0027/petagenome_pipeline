#!/bin/bash

fq1=../../test/s_6_1.fastq.gz
fq2=../../test/s_6_2.fastq.gz

wdir=work
odir=results
refdir=ref

log=t.log

ret=0

wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir} ${wdir}

wfq1=${wdir}/$(basename ${fq1} | sed 's#.gz$##')
wfq2=${wdir}/$(basename ${fq2} | sed 's#.gz$##')

gunzip ${fq1} -c > ${wfq1}
gunzip ${fq2} -c > ${wfq2}

/usr/local/bin/apptainer exec --bind ${wfq1},${wfq2},${odir} ../prinseq.sif prinseq-lite.pl -trim_right 10 -trim_left 10 -trim_qual_right 20 -trim_qual_left 20 -trim_qual_window 20 -min_len 75 -derep 1 -lc_method dust -lc_threshold 7 -trim_ns_right 1 -ns_max_n 0 -fastq ${wfq1} -fastq2 ${wfq2} -out_good ${odir}/good -out_bad ${odir}/bad > ${log} 2>&1
for i in $(ls $odir/*.fastq)
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

# #{prinseq_dir}/bin/prinseq-lite.pl -trim_right 10 -trim_left 10 -trim_qual_right 20 -trim_qual_left 20 -trim_qual_window 20 -min_len 75 -derep 1 -lc_method dust -lc_threshold 7 -trim_ns_right 1 -ns_max_n 0 -fastq #{out_cutadapt}_1.fastq -fastq2 #{out_cutadapt}_2.fastq -out_good #{out_prinseq} -out_bad #{out_prinseq}_bad
# 01_trim_qc.HiSeq.rb




