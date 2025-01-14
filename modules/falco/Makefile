all : falco.sbx falco.sif

clean :
	sudo rm -rf falco.sbx falco.sif build-temp-*

falco.sbx : falco.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox falco.sbx falco.def

falco.sif : falco.sbx
	sudo /usr/local/bin/singularity build --fix-perms falco.sif falco.sbx

shell : falco.sbx
	sudo /usr/local/bin/singularity shell --writable falco.sbx

help : falco.sif
	/usr/local/bin/singularity exec --fakeroot falco.sif falco -h

version : falco.sif
	/usr/local/bin/singularity exec --fakeroot falco.sif falco -v
