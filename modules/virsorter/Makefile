all : virsorter.sbx virsorter.sif

clean :
	sudo rm -rf virsorter.sbx virsorter.sif build-temp-*

virsorter.sbx : virsorter.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox virsorter.sbx virsorter.def

virsorter.sif : virsorter.sbx
	sudo /usr/local/bin/singularity build --fix-perms virsorter.sif virsorter.sbx

shell : virsorter.sbx
	sudo /usr/local/bin/singularity shell --writable virsorter.sbx

hello : virsorter.sif
	/usr/local/bin/singularity exec --fakeroot virsorter.sif echo hello

