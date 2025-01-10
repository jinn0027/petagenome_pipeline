all : cutadapt.sbx cutadapt.sif

clean :
	sudo rm -rf cutadapt.sbx cutadapt.sif build-temp-*

cutadapt.sbx : cutadapt.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox cutadapt.sbx cutadapt.def
	#singularity build --fakeroot --sandbox cutadapt.sbx cutadapt.def

cutadapt.sif : cutadapt.sbx
	sudo /usr/local/bin/singularity build --fix-perms cutadapt.sif cutadapt.sbx
	#singularity build --fakeroot cutadapt.sif cutadapt.sbx

test1 : cutadapt.sbx
	singularity shell cutadapt.sbx

test2 : cutadapt.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work cutadapt.sif

test3 : cutadapt.sif
	/usr/local/bin/singularity -d shell --fakeroot cutadapt.sif

test4 : cutadapt.sif
	/usr/local/bin/singularity exec --fakeroot cutadapt.sif echo hello

