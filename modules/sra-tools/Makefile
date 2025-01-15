all : sra-tools.sbx sra-tools.sif

clean :
	sudo rm -rf sra-tools.sbx sra-tools.sif build-temp-*

sra-tools.sbx : sra-tools.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox sra-tools.sbx sra-tools.def

sra-tools.sif : sra-tools.sbx
	sudo /usr/local/bin/singularity build --fix-perms sra-tools.sif sra-tools.sbx

shell : sra-tools.sbx
	sudo /usr/local/bin/singularity shell --writable sra-tools.sbx

help : sra-tools.sif
	/usr/local/bin/singularity exec --fakeroot sra-tools.sif sratools -h

version : sra-tools.sif
	/usr/local/bin/singularity exec --fakeroot sra-tools.sif sratools -h
