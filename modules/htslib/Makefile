all : htslib.sbx htslib.sif

clean :
	sudo rm -rf htslib.sbx htslib.sif build-temp-*

htslib.sbx : htslib.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox htslib.sbx htslib.def

htslib.sif : htslib.sbx
	sudo /usr/local/bin/singularity build --fix-perms htslib.sif htslib.sbx

shell : htslib.sbx
	sudo /usr/local/bin/singularity shell --writable htslib.sbx

hello : htslib.sif
	/usr/local/bin/singularity exec --fakeroot htslib.sif echo hello

