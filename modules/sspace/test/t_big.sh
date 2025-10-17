#!/bin/bash

n_threads=128


for eid in 1620255 1620256 1620257 1620258
do
    bam=../../samtools/test/results.${eid}/out.bam
    contig=../../megahit/test/results.${eid}/out/final.contigs.fa

    wdir=work
    odir=results.${eid}
    log=t.log.${eid}

    mkdir -p ${wdir} ${odir}
    rm -rf ${wdir}/* ${odir}/* ${log}

    bam=$(cd $(dirname ${bam}) && pwd)/$(basename ${bam})
    wdir=$(cd $(dirname ${wdir}) && pwd)/$(basename ${wdir})
    odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})
    contig=$(cd $(dirname ${contig}) && pwd)/$(basename ${contig})

    # lib_name read_type file_1(BAM) file_2 insert_size(avg_ins) std_dev
    #echo "lib1 ${bam} . 500 50" > ${wdir}/ss.lib
    echo "lib1 /dev/shm/petagenome_pipeline/modules/samtools/test/results.1620255/out . 500 50 FR" > ${wdir}/ss.lib

    date > ${log}
    apptainer exec --bind ${contig},${bam},${wdir} \
        ../sspace.sif sh -c "\
            SSPACE_Basic.pl \
            -l ${wdir}/ss.lib \
            -s ${contig} \
            -k 5 \
            -x 1 \
            -T ${n_threads} \
            -o 5 \
            -z ; \
        " >> ${log} 2>&1
    date >> ${log}
    break
done

exit


