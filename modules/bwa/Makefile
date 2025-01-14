all : bwa.sbx bwa.sif

clean :
	sudo rm -rf bwa.sbx bwa.sif build-temp-*

bwa.sbx : bwa.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox bwa.sbx bwa.def
	#singularity build --fakeroot --sandbox bwa.sbx bwa.def

bwa.sif : bwa.sbx
	sudo /usr/local/bin/singularity build --fix-perms bwa.sif bwa.sbx
	#singularity build --fakeroot bwa.sif bwa.sbx

test1 : bwa.sbx
	singularity shell bwa.sbx

test2 : bwa.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work bwa.sif

test3 : bwa.sif
	/usr/local/bin/singularity -d shell --fakeroot bwa.sif

test4 : bwa.sif
	/usr/local/bin/singularity exec --fakeroot bwa.sif echo hello

