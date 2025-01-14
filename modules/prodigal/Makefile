all : prodigal.sbx prodigal.sif

clean :
	sudo rm -rf prodigal.sbx prodigal.sif build-temp-*

prodigal.sbx : prodigal.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox prodigal.sbx prodigal.def
	#singularity build --fakeroot --sandbox prodigal.sbx prodigal.def

prodigal.sif : prodigal.sbx
	sudo /usr/local/bin/singularity build --fix-perms prodigal.sif prodigal.sbx
	#singularity build --fakeroot prodigal.sif prodigal.sbx

test1 : prodigal.sbx
	singularity shell prodigal.sbx

test2 : prodigal.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work prodigal.sif

test3 : prodigal.sif
	/usr/local/bin/singularity -d shell --fakeroot prodigal.sif

test4 : prodigal.sif
	/usr/local/bin/singularity exec --fakeroot prodigal.sif echo hello

