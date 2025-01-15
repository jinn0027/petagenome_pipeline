all : diamond.sbx diamond.sif

clean :
	sudo rm -rf diamond.sbx diamond.sif build-temp-*

diamond.sbx : diamond.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox diamond.sbx diamond.def

diamond.sif : diamond.sbx
	sudo /usr/local/bin/singularity build --fix-perms diamond.sif diamond.sbx

shell : diamond.sbx
	sudo /usr/local/bin/singularity shell --writable diamond.sbx

help : diamond.sif
	/usr/local/bin/singularity exec --fakeroot diamond.sif diamond help

version : diamond.sif
	/usr/local/bin/singularity exec --fakeroot diamond.sif diamond version

