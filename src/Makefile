# SPPARKS multiple-machine Makefile

SHELL = /bin/sh
#.IGNORE:

# Definitions

ROOT =	spk
EXE =	$(ROOT)_$@
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# Targets

help:
	@echo 'Type "make target" where target is one of:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	for file in $$files; do head -1 $$file; done
	@echo ''
	@echo 'or else:'
	@echo ''
	@echo 'make help = usage information'
	@echo 'make clean = remove all Obj directories and *.o files'
	@echo 'make makelist = create Makefile.list with explicit file list'
	@echo 'make makelib = create Makefile.lib to build library'

clean:
	rm -r Obj_*

tar:
	mkdir spk_src
	rm -f *~ MAKE/*~ STUBS/*~ STUBS/*.o STUBS/libmpi.a
	cp -r *.cpp *.h Makefile* MAKE STUBS in.test in.abc log_abc.spk spk_src
	tar -zcvf spk_src.tgz spk_src
	rm -rf spk_src	

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp -p $(SRC) $(INC) Obj_$@
	@cp MAKE/Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE)  "OBJ = $(OBJ)" "INC = $(INC)" "EXE = ../$(EXE)" ../$(EXE)
	@if [ -d Obj_$@ ]; then cd Obj_$@; rm $(SRC) $(INC) Makefile*; fi

# Update Makefile.lib and Makefile.list

makelib:
	@csh Make.csh Makefile.lib

makelist:
	@csh Make.csh Makefile.list
