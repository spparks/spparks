# g++ = RedHat Linux box, g++, MPICH

SHELL = /bin/sh

# ---------------------------------------------------------------------
# compiler/linker settings
# specify flags and libraries needed for your compiler

CC =		${MPI_HOME}/bin/mpicxx
C =		${MPI_HOME}/bin/mpicc
CCFLAGS =	-g -O -std=c++11
SHFLAGS =	-fPIC
DEPFLAGS =	-M

LINK =		${CC}
LINKFLAGS =	-g -O
LIB =	  	
SIZE =		size

ARCHIVE =	ar
ARFLAGS =	-rc
SHLIBFLAGS =	-shared

# ---------------------------------------------------------------------
# SPPARKS-specific settings
# specify settings for SPPARKS features you will use

# SPPARKS ifdef options, see doc/Section_start.html

SPK_INC =	-DSPPARKS_GZIP  -DSPPARKS_JPEG -DSPPARKS_BIGBIG
SPK_INC =	-DSPPARKS_GZIP  -DSPPARKS_JPEG -DSTITCH_PARALLEL
SPK_INC =	-DSPPARKS_GZIP  -DSPPARKS_JPEG -DSPPARKS_UNORDERED_MAP


# MPI library, can be src/STUBS dummy lib
# INC = path for mpi.h, MPI compiler settings
# PATH = path for MPI library
# LIB = name of MPI library

MPI_INC =       -I${MPI_HOME}/include
MPI_PATH =      -L${MPI_HOME}/lib
MPI_LIB =	-lmpi
MPI_INC = 
MPI_PATH =
MPI_LIB =

# JPEG library, only needed if -DLAMMPS_JPEG listed with LMP_INC
# INC = path for jpeglib.h
# PATH = path for JPEG library
# LIB = name of JPEG library

JPG_INC = -I/usr/local/include
JPG_PATH = -L/usr/local/lib
JPG_LIB = -ljpeg

#STITCH_INC = -I/home/jamitch/local/stitch/include
#STITCH_PATH = -L/home/jamitch/local/stitch/lib
#STITCH_LIB = -lstitch -ldl -lpthread

# ---------------------------------------------------------------------
# build rules and dependencies
# no need to edit this section

include	Makefile.package.settings
include	Makefile.package

EXTRA_INC = $(SPK_INC) $(PKG_INC) $(MPI_INC) $(JPG_INC) $(PKG_SYSINC)
EXTRA_PATH = $(PKG_PATH) $(MPI_PATH) $(JPG_PATH) $(PKG_SYSPATH)
EXTRA_LIB = $(PKG_LIB) $(MPI_LIB) $(JPG_LIB) $(PKG_SYSLIB)
EXTRA_CPP_DEPENDS = $(PKG_CPP_DEPENDS)
EXTRA_LINK_DEPENDS = $(PKG_LINK_DEPENDS)

# Path to src files

vpath %.cpp ..
vpath %.h ..

# Link target

$(EXE):	$(OBJ)
	$(LINK) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(EXTRA_LIB) $(LIB) -o $(EXE)
	$(SIZE) $(EXE)

# Library targets

lib:	$(OBJ)
	$(ARCHIVE) $(ARFLAGS) $(EXE) $(OBJ)

shlib:	$(OBJ)
	$(CC) $(CCFLAGS) $(SHFLAGS) $(SHLIBFLAGS) $(EXTRA_PATH) -o $(EXE) \
        $(OBJ) $(EXTRA_LIB) $(LIB)

# Compilation rules

%.o:%.cpp
	$(CC) $(CCFLAGS) $(SHFLAGS) $(EXTRA_INC) -c $<

%.d:%.cpp
	$(CC) $(CCFLAGS) $(EXTRA_INC) $(DEPFLAGS) $< > $@

# Individual dependencies

depend : fastdep.exe $(SRC)
	@./fastdep.exe $(EXTRA_INC) -- $^ > .depend || exit 1

fastdep.exe: ../DEPEND/fastdep.c
	cc -O -o $@ $<

sinclude .depend
