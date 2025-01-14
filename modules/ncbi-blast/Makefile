all : ncbi-blast.sbx ncbi-blast.sif

clean :
	sudo rm -rf ncbi-blast.sbx ncbi-blast.sif build-temp-*

ncbi-blast.sbx : ncbi-blast.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox ncbi-blast.sbx ncbi-blast.def
	#singularity build --fakeroot --sandbox ncbi-blast.sbx ncbi-blast.def

ncbi-blast.sif : ncbi-blast.sbx
	sudo /usr/local/bin/singularity build --fix-perms ncbi-blast.sif ncbi-blast.sbx
	#singularity build --fakeroot ncbi-blast.sif ncbi-blast.sbx

test1 : ncbi-blast.sbx
	singularity shell ncbi-blast.sbx

test2 : ncbi-blast.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work ncbi-blast.sif

test3 : ncbi-blast.sif
	/usr/local/bin/singularity -d shell --fakeroot ncbi-blast.sif

test4 : ncbi-blast.sif
	/usr/local/bin/singularity exec --fakeroot ncbi-blast.sif echo hello

