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

# List of all targets

help:
	@echo ''
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'
	@echo 'make tar                 spk_src.tar.gz of src dir'
	@echo 'make tar-full            spk_full.tgz of repository trunk'
	@echo 'make makelib             update Makefile.lib for library build'
	@echo 'make makelist            update Makefile.list used by old makes'
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
	@$(SHELL) Make.sh style
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

# Create a tarball of working copy of src dir

tar:
	@cd STUBS; make clean
	@cd ..; tar cvzf src/$(ROOT)_src.tar.gz \
	  src/Make* src/MAKE src/*.cpp src/*.h src/STUBS \
	  --exclude=*/.svn
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
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.lib

makelist:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.list

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

