all : bbmap.sbx bbmap.sif

clean :
	sudo rm -rf bbmap.sbx bbmap.sif build-temp-*

bbmap.sbx : bbmap.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox bbmap.sbx bbmap.def
	#singularity build --fakeroot --sandbox bbmap.sbx bbmap.def

bbmap.sif : bbmap.sbx
	sudo /usr/local/bin/singularity build --fix-perms bbmap.sif bbmap.sbx
	#singularity build --fakeroot bbmap.sif bbmap.sbx

test1 : bbmap.sbx
	singularity shell bbmap.sbx

test2 : bbmap.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work bbmap.sif

test3 : bbmap.sif
	/usr/local/bin/singularity -d shell --fakeroot bbmap.sif

test4 : bbmap.sif
	/usr/local/bin/singularity exec --fakeroot bbmap.sif echo hello

