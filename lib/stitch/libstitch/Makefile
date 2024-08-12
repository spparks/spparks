# 
# megmcca: Aug 6 2024
# ENVIRONMENT variables CC and CXX initialized in-file as in 
# old Makefile
# 
# jamitch: March 19 2024
# ENVIRONMENT variables CC and CXX should initialized
# make -f Makefile.cee stitch.lib
# make -f Makefile.cee stitch_module
#
CC=mpicc
CXX=mpicxx

stitch.lib: stitch.c stitch.h sqlite3.c sqlite3.h
	CC="${CC}" CXX="${CXX}" ${CC} -std=gnu99 -c stitch.c sqlite3.c -DSTITCH_PARALLEL -lstitch -ldl -lpthread
	ar rc libstitch.a stitch.o sqlite3.o

stitch_test: stitch_test.c stitch.lib
	${CC} -o stitch_test stitch_test.c -DSTITCH_PARALLEL -L. -lstitch -ldl -lpthread

## NEW UNTESTED
stitch_module: stitch.c stitch.h sqlite3.c sqlite3.h
	cd .. ; CC="gcc" CXX="g++" python3 setup.py build && CC="gcc" CXX="g++" python3 setup.py install

clean:
	rm -f *.{o,so} libstitch.a stitch_test *.st

