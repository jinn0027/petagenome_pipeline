all : bbmap.sbx bbmap.sif

clean :
	sudo rm -rf bbmap.sbx bbmap.sif build-temp-*

bbmap.sbx : bbmap.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox bbmap.sbx bbmap.def

bbmap.sif : bbmap.sbx
	sudo /usr/local/bin/singularity build --fix-perms bbmap.sif bbmap.sbx

shell : bbmap.sbx
	sudo /usr/local/bin/singularity shell --writable bbmap.sbx

hello : bbmap.sif
	/usr/local/bin/singularity exec --fakeroot bbmap.sif echo hello

