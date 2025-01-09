all : prinseq.sbx prinseq.sif

clean :
	sudo rm -rf prinseq.sbx prinseq.sif build-temp-*

prinseq.sbx : prinseq.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox prinseq.sbx prinseq.def
	#singularity build --fakeroot --sandbox prinseq.sbx prinseq.def

prinseq.sif : prinseq.sbx
	sudo /usr/local/bin/singularity build --fix-perms prinseq.sif prinseq.sbx
	#singularity build --fakeroot prinseq.sif prinseq.sbx

test1 : prinseq.sif
	singularity shell prinseq.sif

test2 : prinseq.sif
	@mkdir -p /tmp/work -m 777
	sudo /usr/local/bin/singularity shell --bind /scratch/local,../../external:/external,/tmp/work:/work prinseq.sif

