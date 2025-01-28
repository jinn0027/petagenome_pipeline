#!/bin/bash -eu

fa1=../../test/hoge.fasta
kv1=../../test/rename.tsv

odir=results
refdir=ref

log=t.log

ret=0

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
kv1=$(cd $(dirname ${kv1}) && pwd)/$(basename ${kv1})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
rm -rf ${odir}/*

/usr/local/bin/apptainer exec --bind ${fa1},${kv1},${odir} ../seqkit.sif sh -c "\
    seqkit replace -p '^(\S+)' -r '{kv}' -k ${kv1} ${fa1} > ${odir}/out.replace.fasta" \
> ${log} 2>&1

/usr/local/bin/apptainer exec --bind ${fa1},${odir} ../seqkit.sif sh -c "\
    seqkit grep -n -r -p '\S+01' -o ${odir}/out.grep.fasta ${fa1}" \
>> ${log} 2>&1

failed=""
for i in $(ls ${refdir}/*.fasta)
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

# #{dir_tools_}/bin/seqkit replace --pattern "^(\w+)(gene\d+_gene_\d+-)" --replacement "{kv}:${2}" --kv-file <(awk '{OFS="\t"} {print $2,$1}' #{query_id_table_}) > #{dir_seq_}/VIRSorter_prophages_cat-4.fa
# #{seqkit_dir}/bin/seqkit grep -n --pattern-file <(cut -f 1 #{out_count}.txt) -o #{contig_proph}.fa
# 06_prophage_detection.02_after_processing.rb
