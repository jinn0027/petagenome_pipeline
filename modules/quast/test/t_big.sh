#!/bin/bash

n_threads=$(nproc)
mem=128
mem_in_byte=${mem}000000000

for eid in 1620255 1620256 1620257 1620258
do
    for k in megahit.contig spades.contig spades.scaffold soapdenovo2.scaffold
    do
        odir=results.${eid}.${k}
        log=t.log.${eid}.${k}

        mkdir -p ${odir}
        rm -rf ${odir}/* ${log}

        if [ "${k}" = "megahit.contig" ] ; then
            fa=../../megahit/test/results.${eid}/out/final.contigs.fa
        elif [ "${k}" = "spades.contig" ] ; then
            fa=../../spades/test/results.${eid}.all/contigs.fasta
        elif [ "${k}" = "spades.scaffold" ] ; then
            fa=../../spades/test/results.${eid}.all/scaffolds.fasta
        elif [ "${k}" = "soapdenovo2.scaffold" ] ; then
            fa=../../soapdenovo2/test/results.${eid}/out.scafSeq
        else
            echo "${k} is not supported"
            exit 1
        fi

        if [ ! -f ${fa} ] ; then
            echo "${fa} is not found"
            continue
        fi

        fa=$(cd $(dirname ${fa}) && pwd)/$(basename ${fa})
        odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

        date > ${log}
        apptainer exec --bind ${fa},${odir} ../quast.sif \
                  quast.py -o ${odir} ${fa} >> ${log} 2>&1
        date >> ${log}
    done
done

