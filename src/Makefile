# SPPARKS multiple-machine Makefile

SHELL = /bin/sh
#.IGNORE:

# Definitions

ROOT =	spk
EXE =	$(ROOT)_$@
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# Package variables

PACKAGE =
PACKAGEUC = $(shell perl -e 'printf("%s", uc("$(PACKAGE)"));')
YESDIR = $(shell perl -e 'printf("%s", uc("$(@:yes-%=%)"));')
NODIR  = $(shell perl -e 'printf("%s", uc("$(@:no-%=%)"));')

# Targets

# List of all targets

help:
	@echo ''
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'
	@echo 'make tar                 spk_src.tar.gz of src dir and packages'
	@echo 'make tar-full            spk_full.tgz of repository trunk'
	@echo 'make makelib             update Makefile.lib for library build'
	@echo 'make makelist            update Makefile.list used by old makes'
	@echo ''
	@echo 'make package             list available packages'
	@echo 'make package-status      status of all packages'
	@echo 'make yes-package         install a single package in src dir'
	@echo 'make no-package          remove a single package from src dir'
	@echo 'make yes-all             install all packages in src dir'
	@echo 'make no-all              remove all packages from src dir'
	@echo ''
	@echo 'make package-update      replace src files with package files'
	@echo 'make package-overwrite   replace package files with src files'
	@echo ''
	@echo 'make machine             build SPPARKS where machine is one of:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''

# Build the code

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@csh Make.csh style
	@cp -p *.cpp *.h Obj_$@
	@cp MAKE/Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = $(INC)" "EXE = ../$(EXE)" ../$(EXE)
	@if [ -d Obj_$@ ]; then cd Obj_$@; rm *.cpp *.h Makefile*; fi

# Remove machine-specific object files

clean:
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'

clean-all:
	rm -rf Obj_*

clean-%:
	rm -rf Obj_$(@:clean-%=%)

# Create a tarball of working copy of src dir and packages

tar:
	@cd STUBS; make clean
	@cd ..; tar cvzf src/$(ROOT)_src.tar.gz \
	  src/Make* src/Package.csh src/MAKE src/*.cpp src/*.h src/STUBS \
	  $(patsubst %,src/%,$(PACKAGEUC)) --exclude=*/.svn
	@cd STUBS; make
	@echo "Created $(ROOT)_src.tar.gz"

# Create a tarball of the latest repository copy

tar-full:
	@mkdir tmp_tardir; cd tmp_tardir; \
	svn export svn+ssh://development.sandia.gov/usr/local/svn/spparks; \
	tar -zcvf ../spk_full.tgz spparks/trunk; \
	cd ..; rm -rf tmp_tardir

# Update Makefile.lib and Makefile.list

makelib:
	@csh Make.csh Makefile.lib

makelist:
	@csh Make.csh Makefile.list

# Package management

package:
	@echo 'Standard packages:' $(PACKAGE)
	@echo ''
	@echo 'make package             list available packages'
	@echo 'make package-status      status of all packages'
	@echo 'make yes-package         install a single package in src dir'
	@echo 'make no-package          remove a single package from src dir'
	@echo 'make yes-all             install all packages in src dir'
	@echo 'make no-all              remove all packages from src dir'
	@echo ''
	@echo 'make package-update      replace src files with package files'
	@echo 'make package-overwrite   replace package files with src files'

yes-all:
	@for p in $(PACKAGE); do $(MAKE) yes-$$p; done

no-all:
	@for p in $(PACKAGE); do $(MAKE) no-$$p; done

yes-%:
	@if [ ! -e $(YESDIR) ]; then \
	  echo "Package $(@:yes-%=%) does not exist"; \
	else \
	  echo "Installing package $(@:yes-%=%)"; \
	  cd $(YESDIR); csh -f Install.csh 1; \
	fi;

no-%:
	@if [ ! -e $(NODIR) ]; then \
	  echo "Package $(@:no-%=%) does not exist"; \
	else \
	  echo "Uninstalling package $(@:no-%=%), ignore errors"; \
	  cd $(NODIR); csh -f Install.csh 0; cd ..; $(MAKE) clean-all; \
        fi;

# status = list differences between src and package files
# update = replace src files with newer package files
# overwrite = overwrite package files with newer src files

package-status:
	@for p in $(PACKAGEUC); do csh -f Package.csh $$p status; done

package-update:
	@for p in $(PACKAGEUC); do csh -f Package.csh $$p update; done

package-overwrite:
	@for p in $(PACKAGEUC); do csh -f Package.csh $$p overwrite; done

# The test feature is not documented on purpose
# Run the tests in directory ../test using Makefile.$(TESTMACHINE)
TESTMACHINE = serial

test: 
	@make $(TESTMACHINE)
	@cp MAKE/Makefile.$(TESTMACHINE) ../test/Makefile
	@cd ../test; rm -f test.out; \
	files="in.ising in.membrane in.potts"; \
	for file in $$files; do echo >> test.out;echo >> test.out; \
	echo "*** Testing file $$file ***" >> test.out;echo >> test.out; \
	$(MAKE) "INPUTFILE = $$file" "EXE = ../src/$(ROOT)_$(TESTMACHINE)" test >> test.out 2>&1; \
	done; rm Makefile

