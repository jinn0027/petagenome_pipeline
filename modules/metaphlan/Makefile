all : metaphlan.sbx metaphlan.sif

clean :
	sudo rm -rf metaphlan.sbx metaphlan.sif build-temp-*

metaphlan.sbx : metaphlan.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox metaphlan.sbx metaphlan.def
	#singularity build --fakeroot --sandbox metaphlan.sbx metaphlan.def

metaphlan.sif : metaphlan.sbx
	sudo /usr/local/bin/singularity build --fix-perms metaphlan.sif metaphlan.sbx
	#singularity build --fakeroot metaphlan.sif metaphlan.sbx

test1 : metaphlan.sbx
	singularity shell metaphlan.sbx

test2 : metaphlan.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work metaphlan.sif

test3 : metaphlan.sif
	/usr/local/bin/singularity -d shell --fakeroot metaphlan.sif

test4 : metaphlan.sif
	/usr/local/bin/singularity exec --fakeroot metaphlan.sif echo hello

