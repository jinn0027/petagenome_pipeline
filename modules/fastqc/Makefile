all : fastqc.sif

clean :
	sudo rm -rf fastqc.sbx fastqc.sif build-temp-*

fastqc.sif : fastqc.def
	sudo /usr/local/bin/singularity build fastqc.sif fastqc.def

test1 : fastqc.sif
	singularity shell fastqc.sif

