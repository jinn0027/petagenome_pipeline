all : prinseq.sbx prinseq.sif

clean :
	sudo rm -rf prinseq.sbx prinseq.sif build-temp-*

prinseq.sbx : prinseq.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox prinseq.sbx prinseq.def

prinseq.sif : prinseq.sbx
	sudo /usr/local/bin/singularity build --fix-perms prinseq.sif prinseq.sbx

shell : prinseq.sbx
	sudo /usr/local/bin/singularity shell --writable prinseq.sbx

help : prinseq.sif
	/usr/local/bin/singularity exec --fakeroot prinseq.sif prinseq-lite.pl -help

version : prinseq.sif
	/usr/local/bin/singularity exec --fakeroot prinseq.sif prinseq-lite.pl -version

