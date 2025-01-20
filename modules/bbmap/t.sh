#!/bin/bash

fq1=../test/s_6_1.fastq.gz
fq2=../test/s_6_2.fastq.gz
odir=results

fq1=$(cd $(dirname ${fq1}) && pwd)/$(basename ${fq1})
fq2=$(cd $(dirname ${fq2}) && pwd)/$(basename ${fq2})
odir=$(cd $(dirname ${odir}) && pwd)/$(basename ${odir})

mkdir -p ${odir}
/usr/local/bin/apptainer exec --bind ${fq1},${fq2},${odir} bbmap.sif bbmap -h

exit

# #{bbmap_dir}/bbmap/bbmap.sh -Xmx#{mem}g threads=#{n_threads} ref=#{contig_}
# #{bbmap_dir}/bbmap/bbmap.sh -Xmx#{mem}g threads=#{n_threads} in=#{query1_} in2=#{query2_} ambiguous=random minid=0.95 pairlen=1500 out=#{out}.sam scafstats=#{out}.scafstats
#10_contig_coverage.rb













