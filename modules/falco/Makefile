all : falco.sbx falco.sif

clean :
	sudo rm -rf falco.sbx falco.sif build-temp-*

falco.sbx : falco.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox falco.sbx falco.def
	#singularity build --fakeroot --sandbox falco.sbx falco.def

falco.sif : falco.sbx
	sudo /usr/local/bin/singularity build --fix-perms falco.sif falco.sbx
	#singularity build --fakeroot falco.sif falco.sbx

test1 : falco.sbx
	singularity shell falco.sbx

test2 : falco.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work falco.sif

test3 : falco.sif
	/usr/local/bin/singularity -d shell --fakeroot falco.sif

test4 : falco.sif
	/usr/local/bin/singularity exec --fakeroot falco.sif echo hello

