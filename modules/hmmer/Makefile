all : hmmer.sbx hmmer.sif

clean :
	sudo rm -rf hmmer.sbx hmmer.sif build-temp-*

hmmer.sbx : hmmer.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox hmmer.sbx hmmer.def

hmmer.sif : hmmer.sbx
	sudo /usr/local/bin/singularity build --fix-perms hmmer.sif hmmer.sbx

shell : hmmer.sbx
	sudo /usr/local/bin/singularity shell --writable hmmer.sbx

help : hmmer.sif
	/usr/local/bin/singularity exec --fakeroot hmmer.sif jackhmmer -h

version : hmmer.sif
	/usr/local/bin/singularity exec --fakeroot hmmer.sif jackhmmer -h | head -n 2

