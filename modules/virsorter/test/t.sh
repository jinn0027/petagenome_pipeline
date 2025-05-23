#!/bin/bash

# Since #cpus may affect the result, it should be fixed here.
n_cpus=1

fa1=../../test/1seq.fa
locdir=/opt/VirSorter
extdir=../../../external

db=/opt/VirSorter/virsorter-data

odir=results
refdir=ref

log=t.log

ret=0

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})
extdir=$(cd $(dirname ${extdir}) && pwd)/$(basename ${extdir})

odir_blst=${odir}/blast
odir_dmnd=${odir}/diamond

mkdir -p ${odir_blst} ${odir_dmnd}
rm -rf ${odir_blst}/* ${odir_dmnd}/*

apptainer exec --bind ${fa1},${odir_blst},${extdir}/virsorter-data:${locdir}/virsorter-data,${extdir}/mga_linux_ia64:${locdir}/mga_linux_ia64 ../virsorter.sbx wrapper_phage_contigs_sorter_iPlant.pl --db 1 --data-dir ${db} --ncpu ${n_cpus} --wdir ${odir_blst} --fna ${fa1} > ${log} 2>&1

apptainer exec --bind ${fa1},${odir_dmnd},${extdir}/virsorter-data:${locdir}/virsorter-data,${extdir}/mga_linux_ia64:${locdir}/mga_linux_ia64 ../virsorter.sbx wrapper_phage_contigs_sorter_iPlant.pl --diamond --db 2 -data-dir ${db} -ncpu ${n_cpus} -wdir ${odir_dmnd} --fna ${fa1} >> ${log} 2>&1

failed=""
for i in $(ls ${refdir}/blast/*.csv)
do
    j=${odir_blst}/$(basename $i)
    diff -q $i $j >> ${log} 2>&1 && :
    if [ $? -ne 0 ]; then
	failed="${failed} blast/$(basename $i)"
        ret=1
    fi
done
for i in $(ls ${refdir}/diamond/*.csv)
do
    j=${odir_dmnd}/$(basename $i)
    diff -q $i $j >> ${log} 2>&1 && :
    if [ $? -ne 0 ]; then
	failed="${failed} diamond/$(basename $i)"
        ret=1
    fi
done

if [ ${ret} -eq 0 ]; then
    echo " PASSED" | tee -a ${log}
else
    echo " FAILED : ${failed}" | tee -a ${log}
fi

exit ${ret}

# perl #{virsorter_dir}/wrapper_phage_contigs_sorter_iPlant.pl #{opt_virome} --db #{vs_db_id} -data-dir #{dir_data_vs_} -ncpu #{n_threads} -wdir #{dir_out_} -f #{query_fa_}
# 06_prophage_detection.01_each_sample.rb












