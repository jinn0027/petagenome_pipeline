all : bedtools.sbx bedtools.sif

clean :
	sudo rm -rf bedtools.sbx bedtools.sif build-temp-*

bedtools.sbx : bedtools.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox bedtools.sbx bedtools.def
	#singularity build --fakeroot --sandbox bedtools.sbx bedtools.def

bedtools.sif : bedtools.sbx
	sudo /usr/local/bin/singularity build --fix-perms bedtools.sif bedtools.sbx
	#singularity build --fakeroot bedtools.sif bedtools.sbx

test1 : bedtools.sbx
	singularity shell bedtools.sbx

test2 : bedtools.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work bedtools.sif

test3 : bedtools.sif
	/usr/local/bin/singularity -d shell --fakeroot bedtools.sif

test4 : bedtools.sif
	/usr/local/bin/singularity exec --fakeroot bedtools.sif echo hello

