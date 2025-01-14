all : prodigal.sbx prodigal.sif

clean :
	sudo rm -rf prodigal.sbx prodigal.sif build-temp-*

prodigal.sbx : prodigal.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox prodigal.sbx prodigal.def

prodigal.sif : prodigal.sbx
	sudo /usr/local/bin/singularity build --fix-perms prodigal.sif prodigal.sbx

shell : prodigal.sbx
	sudo /usr/local/bin/singularity shell --writable prodigal.sbx

help : prodigal.sif
	/usr/local/bin/singularity exec --fakeroot prodigal.sif prodigal -h

version : prodigal.sif
	/usr/local/bin/singularity exec --fakeroot prodigal.sif prodigal -v



