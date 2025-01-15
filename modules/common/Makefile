all : el9.sbx el9.sif

clean :
	sudo rm -rf el9.sbx el9.sif build-temp-*

el9.sbx : el9.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox el9.sbx el9.def

el9.sif : el9.sbx
	sudo /usr/local/bin/singularity build --fix-perms el9.sif el9.sbx

shell : el9.sbx
	sudo /usr/local/bin/singularity shell --writable el9.sbx

help : el9.sif
	/usr/local/bin/singularity exec --fakeroot el9.sif fastqc --help

version : el9.sif
	/usr/local/bin/singularity exec --fakeroot el9.sif fastqc --version

