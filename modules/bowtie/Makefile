all : bowtie.sbx bowtie.sif

clean :
	sudo rm -rf bowtie.sbx bowtie.sif build-temp-*

bowtie.sbx : bowtie.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox bowtie.sbx bowtie.def

bowtie.sif : bowtie.sbx
	sudo /usr/local/bin/singularity build --fix-perms bowtie.sif bowtie.sbx

shell : bowtie.sbx
	sudo /usr/local/bin/singularity shell --writable bowtie.sbx

help : bowtie.sif
	/usr/local/bin/singularity exec --fakeroot bowtie.sif bowtie --help

version : bowtie.sif
	/usr/local/bin/singularity exec --fakeroot bowtie.sif bowtie --version

