#!/bin/bash

fa1=../../test/test.iupac.fa
bed1=../../test/test.iupac.bed
bed2=../../test/precisionTest2.bed
bed3=../../test/x.bed
bed4=../../test/y.bed

wdir=work
odir=results
refdir=ref

log=t.log

ret=0

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
bed1=$(cd $(dirname ${bed1}) && pwd)/$(basename ${bed1})
bed2=$(cd $(dirname ${bed2}) && pwd)/$(basename ${bed2})
bed3=$(cd $(dirname ${bed3}) && pwd)/$(basename ${bed3})
bed4=$(cd $(dirname ${bed4}) && pwd)/$(basename ${bed4})

wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir} ${wdir}
/usr/local/bin/apptainer exec --bind ${fa1},${bed1},${wdir},${odir} ../bedtools2.sif sh -c "bedtools getfasta -fi ${fa1} -bed ${bed1} -fo ${odir}/picked.fa" > ${log} 2>&1

/usr/local/bin/apptainer exec --bind ${bed2},${wdir},${odir} ../bedtools2.sif sh -c "bedtools sort -i ${bed2} > ${wdir}/sorted.bed && bedtools merge -i ${wdir}/sorted.bed -d 200000 -s -c 4,5,6,7,8 -o distinct,distinct,distinct,distinct,min > ${odir}/merged.bed" >> ${log} 2>&1

/usr/local/bin/apptainer exec --bind ${bed3},${bed4},${wdir},${odir} ../bedtools2.sif sh -c "bedtools intersect -a ${bed3} -b ${bed4} > ${odir}/intersected.bed" >> ${log} 2>&1

for i in $(ls ${refdir}/*.fa ${refdir}/*.bed)
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

# #{bedtools_dir}/bin/bedtools merge -i #{blast_out}.bed -d 150 -s -c 4,5,6,7,8,9 -o distinct,distinct,distinct,distinct,min,max
# #{bedtools_dir}/bin/bedtools merge -i #{blast_out}.merged.tmp2.bed -d 0 -c 4,5,6 -o distinct,distinct,distinct
# #{bedtools_dir}/bin/bedtools getfasta -fi #{contig_proph}.fa -bed #{blast_out}.merged.bed -fo #{merged_prophage}.original.fa
# #{bedtools_dir}/bin/bedtools intersect -a #{orf_position}.bed -b #{read_cov}.region.bed -wo > #{out_orf}.intersect.bed
# 06_prophage_detection.02_after_processing.rb 10_orf_coverage.rb











