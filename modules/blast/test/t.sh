#!/bin/bash

n_threads=$(nproc)
n_al="1"
pi_thre="95"
e_thre="1e-20"

fa1=../../test/ecoli_1K_1.fa.gz
fa2=../../test/q.fa

wdir=work
odir=results
refdir=ref

log=t.log

ret=0

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
fa2=$(cd $(dirname ${fa2}) && pwd)/$(basename ${fa2})
wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir} ${wdir}
rm -rf ${odir}/* ${wdir}/*

wfa1=${wdir}/$(basename ${fa1} | sed 's#.gz$##')
db=${wdir}/$(basename ${fa1})

gunzip ${fa1} -c > ${wfa1}

apptainer exec --bind ${wdir} ../blast.sbx makeblastdb -in ${wfa1} -out ${db} -dbtype nucl -parse_seqids > ${log} 2>&1

apptainer exec --bind ${wdir},${fa2},${odir} ../blast.sbx blastn -task megablast -num_threads ${n_threads} -query ${fa2} -db ${db} -perc_identity ${pi_thre} -evalue ${e_thre} -outfmt 6 -num_alignments ${n_al} -out ${odir}/out.txt >> ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.txt)
do
    j=${odir}/$(basename $i)
    diff -q $i $j >> ${log} 2>&1 && :
    if [ $? -ne 0 ]; then
        failed="${failed} $(basename $i)"
        ret=1
    fi
done

if [ ${ret} -eq 0 ]; then
    echo " PASSED" | tee -a ${log}
else
    echo " FAILED : ${failed}" | tee -a ${log}
fi

exit ${ret}

# #{blast_dir}/makeblastdb -in #{out}.1000.fa -out #{out}.1000 -dbtype nucl -parse_seqids
# #{blast_dir}/bin/blastn -task megablast -num_threads #{n_threads} -query #{prophage}.fa -db #{contig_proph} -perc_identity #{pi_thre} -evalue #{e_thre} -outfmt 6 -num_alignments #{n_al} -out #{blast_out}.tx
# 03_assembly.rb 04_pool_contigs.rb 05_circular_contigs.rb 06_prophage_detection.02_after_processing.rb











