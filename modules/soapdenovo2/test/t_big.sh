#!/bin/bash

n_threads=128
k_mer=41

for eid in 1620255 1620256 1620257 1620258
do
    fq1=/scratch/local/data/metagenome/ERR${eid}_XXXXXXXX_XXXXXXXX_L001_R1_001.fastq.gz
    fq2=/scratch/local/data/metagenome/ERR${eid}_XXXXXXXX_XXXXXXXX_L001_R2_001.fastq.gz
    contig=../../megahit/test/results.${eid}/out/final.contigs.fa
    config_templ=./config.templ

    odir=results.${eid}
    log=t.log.${eid}

    fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
    fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
    contig=$(cd $(dirname ${contig}) && pwd)/$(basename ${contig})
    config_templ=$(cd $(dirname ${config_templ}) && pwd)/$(basename ${config_templ})

    odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

    mkdir -p ${odir}
    rm -rf ${odir}/* ${log}

    date > ${log}
    apptainer exec --bind ${fq1},${fq2},${contig},${config_templ},${odir} ../soapdenovo2.sif sh -c "\
            cat ${config_templ} > ${odir}/config ; \
            echo \"q1=${fq1}\" >> ${odir}/config ; \
            echo \"q2=${fq2}\" >> ${odir}/config ; \
            SOAPdenovo-fusion -D -s ${odir}/config -p ${n_threads} -K ${k_mer} -g ${odir}/out -c ${contig}; \
            SOAPdenovo-127mer map -s ${odir}/config -p ${n_threads} -g ${odir}/out; \
            SOAPdenovo-127mer scaff -p ${n_threads} -g ${odir}/out; \
        " >> ${log} 2>&1
    date >> ${log}
done
