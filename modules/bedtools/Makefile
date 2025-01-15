all : bedtools.sbx bedtools.sif

clean :
	sudo rm -rf bedtools.sbx bedtools.sif build-temp-*

bedtools.sbx : bedtools.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox bedtools.sbx bedtools.def

bedtools.sif : bedtools.sbx
	sudo /usr/local/bin/singularity build --fix-perms bedtools.sif bedtools.sbx

shell : bedtools.sbx
	sudo /usr/local/bin/singularity shell --writable bedtools.sbx

help : bedtools.sif
	/usr/local/bin/singularity exec --fakeroot bedtools.sif bedtools --help

version : bedtools.sif
	/usr/local/bin/singularity exec --fakeroot bedtools.sif bedtools --version
