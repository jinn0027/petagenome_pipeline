all : seqkit.sbx seqkit.sif

clean :
	sudo rm -rf seqkit.sbx seqkit.sif build-temp-*

seqkit.sbx : seqkit.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox seqkit.sbx seqkit.def
	#singularity build --fakeroot --sandbox seqkit.sbx seqkit.def

seqkit.sif : seqkit.sbx
	sudo /usr/local/bin/singularity build --fix-perms seqkit.sif seqkit.sbx
	#singularity build --fakeroot seqkit.sif seqkit.sbx

test1 : seqkit.sbx
	singularity shell seqkit.sbx

test2 : seqkit.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work seqkit.sif

test3 : seqkit.sif
	/usr/local/bin/singularity -d shell --fakeroot seqkit.sif

test4 : seqkit.sif
	/usr/local/bin/singularity exec --fakeroot seqkit.sif echo hello

