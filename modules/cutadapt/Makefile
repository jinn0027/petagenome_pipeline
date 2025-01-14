all : cutadapt.sbx cutadapt.sif

clean :
	sudo rm -rf cutadapt.sbx cutadapt.sif build-temp-*

cutadapt.sbx : cutadapt.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox cutadapt.sbx cutadapt.def

cutadapt.sif : cutadapt.sbx
	sudo /usr/local/bin/singularity build --fix-perms cutadapt.sif cutadapt.sbx

shell : cutadapt.sbx
	sudo /usr/local/bin/singularity shell --writable cutadapt.sbx

help : cutadapt.sif
	/usr/local/bin/singularity exec --fakeroot cutadapt.sif cutadapt --help

version : cutadapt.sif
	/usr/local/bin/singularity exec --fakeroot cutadapt.sif cutadapt --version
