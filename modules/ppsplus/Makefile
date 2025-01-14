all : ppsplus.sbx ppsplus.sif

clean :
	sudo rm -rf ppsplus.sbx ppsplus.sif build-temp-*

ppsplus.sbx : ppsplus.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox ppsplus.sbx ppsplus.def
	#singularity build --fakeroot --sandbox ppsplus.sbx ppsplus.def

ppsplus.sif : ppsplus.sbx
	sudo /usr/local/bin/singularity build --fix-perms ppsplus.sif ppsplus.sbx
	#singularity build --fakeroot ppsplus.sif ppsplus.sbx

test1 : ppsplus.sbx
	singularity shell ppsplus.sbx

test2 : ppsplus.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work ppsplus.sif

test3 : ppsplus.sif
	/usr/local/bin/singularity -d shell --fakeroot ppsplus.sif

test4 : ppsplus.sif
	/usr/local/bin/singularity exec --fakeroot ppsplus.sif echo hello

