# mpicc and mpicxx must be on your path; Set these appropriately

CC=mpicc
CXX=mpicxx

stitch.lib: stitch.c stitch.h sqlite3.c sqlite3.h
	CC="${CC}" CXX="${CXX}" ${CC} -std=gnu99 -c stitch.c sqlite3.c -DSTITCH_PARALLEL -lstitch -ldl -lpthread
	ar rc libstitch.a stitch.o sqlite3.o

stitch_test: stitch_test.c stitch.lib
	${CC} -std=c99 -o stitch_test stitch_test.c -DSTITCH_PARALLEL -L. -lstitch -ldl -lpthread

stitch_module: stitch.c stitch.h sqlite3.c sqlite3.h
	cd .. ; CC="gcc" CXX="g++" python3 setup.py build && CC="gcc" CXX="g++" python3 setup.py install

clean:
	rm -f *.o libstitch.a *.so 
	cd .. ; python3 setup.py clean
	rm -rf ../build
	find .. -name __pycache__ -exec rm -rf {} \; ; find .. -name "*.pyc" -exec rm -rf {} \; ; rm -f *.st
