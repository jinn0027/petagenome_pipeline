all : bwa.sbx bwa.sif

clean :
	sudo rm -rf bwa.sbx bwa.sif build-temp-*

bwa.sbx : bwa.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox bwa.sbx bwa.def

bwa.sif : bwa.sbx
	sudo /usr/local/bin/singularity build --fix-perms bwa.sif bwa.sbx

shell : bwa.sbx
	sudo /usr/local/bin/singularity shell --writable bwa.sbx

help : bwa.sif
	/usr/local/bin/singularity exec --fakeroot bwa.sif bwa

version : bwa.sif
	/usr/local/bin/singularity exec --fakeroot bwa.sif bwa | head -n 2


