all : virsorter.sbx virsorter.sif

clean :
	sudo rm -rf virsorter.sbx virsorter.sif build-temp-*

virsorter.sbx : virsorter.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox virsorter.sbx virsorter.def
	#singularity build --fakeroot --sandbox virsorter.sbx virsorter.def

virsorter.sif : virsorter.sbx
	sudo /usr/local/bin/singularity build --fix-perms virsorter.sif virsorter.sbx
	#singularity build --fakeroot virsorter.sif virsorter.sbx

test1 : virsorter.sbx
	singularity shell virsorter.sbx

test2 : virsorter.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work virsorter.sif

test3 : virsorter.sif
	/usr/local/bin/singularity -d shell --fakeroot virsorter.sif

test4 : virsorter.sif
	/usr/local/bin/singularity exec --fakeroot virsorter.sif echo hello

