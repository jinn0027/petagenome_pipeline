all : cdhit.sbx cdhit.sif

clean :
	sudo rm -rf cdhit.sbx cdhit.sif build-temp-*

cdhit.sbx : cdhit.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox cdhit.sbx cdhit.def

cdhit.sif : cdhit.sbx
	sudo /usr/local/bin/singularity build --fix-perms cdhit.sif cdhit.sbx

shell : cdhit.sbx
	sudo /usr/local/bin/singularity shell --writable cdhit.sbx

help : cdhit.sif
	/usr/local/bin/singularity exec --fakeroot cdhit.sif cd-hit --help

version : cdhit.sif
	/usr/local/bin/singularity exec --fakeroot cdhit.sif cd-hit --help | head -n 1

