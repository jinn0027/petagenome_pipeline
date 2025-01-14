all : seqkit.sbx seqkit.sif

clean :
	sudo rm -rf seqkit.sbx seqkit.sif build-temp-*

seqkit.sbx : seqkit.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox seqkit.sbx seqkit.def

seqkit.sif : seqkit.sbx
	sudo /usr/local/bin/singularity build --fix-perms seqkit.sif seqkit.sbx

shell : seqkit.sbx
	sudo /usr/local/bin/singularity shell --writable seqkit.sbx

hello : seqkit.sif
	/usr/local/bin/singularity exec --fakeroot seqkit.sif echo hello

