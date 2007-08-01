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
	@echo 'make tar = make tarball of src dir'
	@echo 'make package = see what packages are available'

clean:
	rm -r Obj_*

tar:
	mkdir spk_src
	rm -f *~ MAKE/*~ STUBS/*~ STUBS/*.o STUBS/libmpi.a
	cp -r *.cpp *.h Makefile* MAKE STUBS in.* spk_src
	cp -r GPPT spk_src
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

# Packages

package:
	@echo 'Available packages:'
	@echo '  gppt'
	@echo 'make yes-name     to include a package'
	@echo 'make no-name      to exclude a package'
	@echo 'make yes-all      to include all packages'
	@echo 'make no-all       to exclude all packages'

yes-all:
	make yes-gppt

no-all:
	@echo 'Removing files, ignore any rm errors ...'
	@cd GPPT; csh -f Install.csh 0
	@make clean

yes-gppt:
	@cd GPPT; csh -f Install.csh 1
no-gppt:
	@echo 'Removing files, ignore any rm errors ...'
	@cd GPPT; csh -f Install.csh 0
	@make clean

# update src files with package files

package-update:
	@csh -f Package.csh GPPT update

# overwrite package files with src files

package-overwrite:
	@csh -f Package.csh GPPT overwrite

# check differences between src and pacakge files

package-check:
	@csh -f Package.csh GPPT check
