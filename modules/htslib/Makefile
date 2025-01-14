all : htslib.sbx htslib.sif

clean :
	sudo rm -rf htslib.sbx htslib.sif build-temp-*

htslib.sbx : htslib.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox htslib.sbx htslib.def
	#singularity build --fakeroot --sandbox htslib.sbx htslib.def

htslib.sif : htslib.sbx
	sudo /usr/local/bin/singularity build --fix-perms htslib.sif htslib.sbx
	#singularity build --fakeroot htslib.sif htslib.sbx

test1 : htslib.sbx
	singularity shell htslib.sbx

test2 : htslib.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work htslib.sif

test3 : htslib.sif
	/usr/local/bin/singularity -d shell --fakeroot htslib.sif

test4 : htslib.sif
	/usr/local/bin/singularity exec --fakeroot htslib.sif echo hello

