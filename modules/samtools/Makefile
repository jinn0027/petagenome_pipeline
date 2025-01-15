all : samtools.sbx samtools.sif

clean :
	sudo rm -rf samtools.sbx samtools.sif build-temp-*

samtools.sbx : samtools.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox samtools.sbx samtools.def

samtools.sif : samtools.sbx
	sudo /usr/local/bin/singularity build --fix-perms samtools.sif samtools.sbx

shell : samtools.sbx
	sudo /usr/local/bin/singularity shell --writable samtools.sbx

help : samtools.sif
	/usr/local/bin/singularity exec --fakeroot samtools.sif samtools help

version : samtools.sif
	/usr/local/bin/singularity exec --fakeroot samtools.sif samtools version
