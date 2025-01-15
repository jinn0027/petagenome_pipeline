all : metaphlan.sbx metaphlan.sif

clean :
	sudo rm -rf metaphlan.sbx metaphlan.sif build-temp-*

metaphlan.sbx : metaphlan.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox metaphlan.sbx metaphlan.def

metaphlan.sif : metaphlan.sbx
	sudo /usr/local/bin/singularity build --fix-perms metaphlan.sif metaphlan.sbx

shell : metaphlan.sbx
	sudo /usr/local/bin/singularity shell --writable metaphlan.sbx

help : metaphlan.sif
	/usr/local/bin/singularity exec --fakeroot metaphlan.sif metaphlan --help

version : metaphlan.sif
	/usr/local/bin/singularity exec --fakeroot metaphlan.sif metaphlan --version

