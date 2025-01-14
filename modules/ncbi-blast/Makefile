all : ncbi-blast.sbx ncbi-blast.sif

clean :
	sudo rm -rf ncbi-blast.sbx ncbi-blast.sif build-temp-*

ncbi-blast.sbx : ncbi-blast.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox ncbi-blast.sbx ncbi-blast.def

ncbi-blast.sif : ncbi-blast.sbx
	sudo /usr/local/bin/singularity build --fix-perms ncbi-blast.sif ncbi-blast.sbx

shell : ncbi-blast.sbx
	sudo /usr/local/bin/singularity shell --writable ncbi-blast.sbx

help : ncbi-blast.sif
	/usr/local/bin/singularity exec --fakeroot ncbi-blast.sif blastn -help

version : ncbi-blast.sif
	/usr/local/bin/singularity exec --fakeroot ncbi-blast.sif blastn -version

