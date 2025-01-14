all : ppsplus.sbx ppsplus.sif

clean :
	sudo rm -rf ppsplus.sbx ppsplus.sif build-temp-*

ppsplus.sbx : ppsplus.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox ppsplus.sbx ppsplus.def

ppsplus.sif : ppsplus.sbx
	sudo /usr/local/bin/singularity build --fix-perms ppsplus.sif ppsplus.sbx

shell : ppsplus.sbx
	sudo /usr/local/bin/singularity shell --writable ppsplus.sbx

hello : ppsplus.sif
	/usr/local/bin/singularity exec --fakeroot ppsplus.sif echo hello

