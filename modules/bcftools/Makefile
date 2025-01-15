all : bcftools.sbx bcftools.sif

clean :
	sudo rm -rf bcftools.sbx bcftools.sif build-temp-*

bcftools.sbx : bcftools.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox bcftools.sbx bcftools.def

bcftools.sif : bcftools.sbx
	sudo /usr/local/bin/singularity build --fix-perms bcftools.sif bcftools.sbx

shell : bcftools.sbx
	sudo /usr/local/bin/singularity shell --writable bcftools.sbx

help : bcftools.sif
	/usr/local/bin/singularity exec --fakeroot bcftools.sif bcftools --help

version : bcftools.sif
	/usr/local/bin/singularity exec --fakeroot bcftools.sif bcftools --version

