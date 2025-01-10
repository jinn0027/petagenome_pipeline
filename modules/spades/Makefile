all : spades.sbx spades.sif

clean :
	sudo rm -rf spades.sbx spades.sif build-temp-*

spades.sbx : spades.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox spades.sbx spades.def
	#singularity build --fakeroot --sandbox spades.sbx spades.def

spades.sif : spades.sbx
	sudo /usr/local/bin/singularity build --fix-perms spades.sif spades.sbx
	#singularity build --fakeroot spades.sif spades.sbx

test1 : spades.sbx
	singularity shell spades.sbx

test2 : spades.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work spades.sif

test3 : spades.sif
	/usr/local/bin/singularity -d shell --fakeroot spades.sif

test4 : spades.sif
	/usr/local/bin/singularity exec --fakeroot spades.sif echo hello

