all : hmmer.sbx hmmer.sif

clean :
	sudo rm -rf hmmer.sbx hmmer.sif build-temp-*

hmmer.sbx : hmmer.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox hmmer.sbx hmmer.def
	#singularity build --fakeroot --sandbox hmmer.sbx hmmer.def

hmmer.sif : hmmer.sbx
	sudo /usr/local/bin/singularity build --fix-perms hmmer.sif hmmer.sbx
	#singularity build --fakeroot hmmer.sif hmmer.sbx

test1 : hmmer.sbx
	singularity shell hmmer.sbx

test2 : hmmer.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work hmmer.sif

test3 : hmmer.sif
	/usr/local/bin/singularity -d shell --fakeroot hmmer.sif

test4 : hmmer.sif
	/usr/local/bin/singularity exec --fakeroot hmmer.sif echo hello

