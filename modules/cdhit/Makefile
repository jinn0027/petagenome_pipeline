all : cdhit.sbx cdhit.sif

clean :
	sudo rm -rf cdhit.sbx cdhit.sif build-temp-*

cdhit.sbx : cdhit.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox cdhit.sbx cdhit.def
	#singularity build --fakeroot --sandbox cdhit.sbx cdhit.def

cdhit.sif : cdhit.sbx
	sudo /usr/local/bin/singularity build --fix-perms cdhit.sif cdhit.sbx
	#singularity build --fakeroot cdhit.sif cdhit.sbx

test1 : cdhit.sbx
	singularity shell cdhit.sbx

test2 : cdhit.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work cdhit.sif

test3 : cdhit.sif
	/usr/local/bin/singularity -d shell --fakeroot cdhit.sif

test4 : cdhit.sif
	/usr/local/bin/singularity exec --fakeroot cdhit.sif echo hello

