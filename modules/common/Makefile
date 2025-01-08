all : el9.sbx

clean :
	sudo rm -rf el9.sbx el9.sif build-temp-*

#el9.sif : el9.def
#	singularity build --fakeroot el9.sif el9.def

el9.sbx : el9.def
	sudo /usr/local/bin/singularity build --sandbox el9.sbx el9.def
	#singularity build --fakeroot --sandbox el9.sbx el9.def

el9.sif : el9.sbx
	sudo /usr/local/bin/singularity build el9.sif el9.sbx
	#singularity build --fakeroot el9.sif el9.sbx

test1 :
	singularity shell el9.sif

#test2 : el9.sbx
#	sudo /usr/local/bin/singularity shell --writable el9.sbx

test3 : el9.sbx
	@mkdir -p /tmp/work -m 777
	#singularity shell --fakeroot --writable --bind /scratch/local,../../external:/external,/tmp/work:/work el9.sbx
	sudo /usr/local/bin/singularity shell --writable --bind /scratch/local,../../external:/external,/tmp/work:/work el9.sbx

test4 : el9.sbx
	@mkdir -p /tmp/work -m 777
	singularity shell --writable --bind /scratch/local,../../external:/external,/tmp/work:/work el9.sbx

test5 : el9.sbx
	@mkdir -p /tmp/work -m 777
	singularity run --writable --bind /scratch/local,../../external:/external,/tmp/work:/work el9.sbx

test6 : el9.sif
	@mkdir -p /tmp/work -m 777
	singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work el9.sif

test7 : el9.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work el9.sif

test8 : el9.sif
	@mkdir -p /tmp/work -m 777
	singularity run --bind /scratch/local,../../external:/external,/tmp/work:/work el9.sif

test9 : el9.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity run --bind /scratch/local,../../external:/external,/tmp/work:/work el9.sif

test10 : el9.sbx
	@mkdir -p /tmp/work -m 777
	singularity run --bind /scratch/local,../../external:/external,/tmp/work:/work el9.sbx

test11 : el9.sbx
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity run --bind /scratch/local,../../external:/external,/tmp/work:/work el9.sbx

test12 : el9.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work el9.sif

