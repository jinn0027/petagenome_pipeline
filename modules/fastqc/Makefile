all : fastqc.sbx fastqc.sif

clean :
	sudo rm -rf fastqc.sbx fastqc.sif build-temp-*

fastqc.sbx : fastqc.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox fastqc.sbx fastqc.def

fastqc.sif : fastqc.sbx
	sudo /usr/local/bin/singularity build --fix-perms fastqc.sif fastqc.sbx

shell : fastqc.sbx
	sudo /usr/local/bin/singularity shell --writable fastqc.sbx

help : fastqc.sif
	/usr/local/bin/singularity exec --fakeroot fastqc.sif fastqc --help

version : fastqc.sif
	/usr/local/bin/singularity exec --fakeroot fastqc.sif fastqc --version
