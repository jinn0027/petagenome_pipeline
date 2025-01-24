#!/bin/bash

n_threads=$(nproc)

#fa1=../../test/ecoli_1K_1.fa.gz
#fa1=../../test/NC_086348.1.fna
#fa1=../../test/NC_091965.1.fna
fa1=../../test/NC_083851.1.fna
db=/opt/VirSorter/virsorter-data

odir=results
refdir=ref

log=t.log

ret=0

fa1=$(cd $(dirname ${fa1}) && pwd)/$(basename ${fa1})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

odir_blst=${odir}/blast
odir_dmnd=${odir}/diamond

mkdir -p ${odir_blst} ${odir_dmnd}

/usr/local/bin/apptainer exec --bind ${fa1},${odir_blst} ../virsorter.sbx wrapper_phage_contigs_sorter_iPlant.pl --db 1 --data-dir ${db} --ncpu ${n_threads} --wdir ${odir_blst} --fna ${fa1} > ${log} 2>&1

/usr/local/bin/apptainer exec --bind ${fa1},${odir_dmnd} ../virsorter.sbx wrapper_phage_contigs_sorter_iPlant.pl --diamond --db 2 -data-dir ${db} -ncpu ${n_threads} -wdir ${odir_dmnd} --fna ${fa1} >> ${log} 2>&1

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












