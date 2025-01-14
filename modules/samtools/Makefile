all : samtools.sbx samtools.sif

clean :
	sudo rm -rf samtools.sbx samtools.sif build-temp-*

samtools.sbx : samtools.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox samtools.sbx samtools.def
	#singularity build --fakeroot --sandbox samtools.sbx samtools.def

samtools.sif : samtools.sbx
	sudo /usr/local/bin/singularity build --fix-perms samtools.sif samtools.sbx
	#singularity build --fakeroot samtools.sif samtools.sbx

test1 : samtools.sbx
	singularity shell samtools.sbx

test2 : samtools.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work samtools.sif

test3 : samtools.sif
	/usr/local/bin/singularity -d shell --fakeroot samtools.sif

test4 : samtools.sif
	/usr/local/bin/singularity exec --fakeroot samtools.sif echo hello

