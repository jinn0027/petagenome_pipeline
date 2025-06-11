.PHONY : all
all : virsorter.sif

virsorter.sif :
	+make -f Makefile.virsorter virsorter.sif

.PHONY : clean
clean :
	+make -f Makefile.virsorter clean
	+make -f Makefile.virsorter_update_db clean

.PHONY : db_updated
db_updated :
	+make -f Makefile.virsorter_update_db db_updated
