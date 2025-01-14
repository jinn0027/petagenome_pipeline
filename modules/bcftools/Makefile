all : bcftools.sbx bcftools.sif

clean :
	sudo rm -rf bcftools.sbx bcftools.sif build-temp-*

bcftools.sbx : bcftools.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox bcftools.sbx bcftools.def
	#singularity build --fakeroot --sandbox bcftools.sbx bcftools.def

bcftools.sif : bcftools.sbx
	sudo /usr/local/bin/singularity build --fix-perms bcftools.sif bcftools.sbx
	#singularity build --fakeroot bcftools.sif bcftools.sbx

test1 : bcftools.sbx
	singularity shell bcftools.sbx

test2 : bcftools.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work bcftools.sif

test3 : bcftools.sif
	/usr/local/bin/singularity -d shell --fakeroot bcftools.sif

test4 : bcftools.sif
	/usr/local/bin/singularity exec --fakeroot bcftools.sif echo hello

