all : spades.sbx spades.sif

clean :
	sudo rm -rf spades.sbx spades.sif build-temp-*

spades.sbx : spades.def
	sudo /usr/local/bin/singularity build --fix-perms --sandbox spades.sbx spades.def

spades.sif : spades.sbx
	sudo /usr/local/bin/singularity build --fix-perms spades.sif spades.sbx

shell : spades.sbx
	sudo /usr/local/bin/singularity shell --writable spades.sbx

help : spades.sif
	/usr/local/bin/singularity exec --fakeroot spades.sif spades.py -h

version : spades.sif
	/usr/local/bin/singularity exec --fakeroot spades.sif spades.py -v
