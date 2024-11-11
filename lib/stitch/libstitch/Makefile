#
# jamitch: Sept 4 2024
# ENVIRONMENT variables CC and CXX should initialized
#    appropriately
# CC=mpicc
# CXX=mpicxx
# 
# Example usage:
# make stitch_test
#     Test runs on 4-procs. Run the test.
#     $ mpirun -np 4 ./stitch_test
#     Expected OUTPUT Below
#     Verifies sanity of build
#     nprocs in x, y, z: 2, 2, 1
#     block size in x, y, z: 2, 2, 3
#     increment in x, y, z: 3, 3, 1
#     total domain bounds (min) - (max): (0, 0, 0) - (13, 13, 4)
#     sweeping 0->n in X, then Y, then Z
#     tests passed
#
# make stitch.lib
# make make clean


stitch.lib: stitch.c stitch.h sqlite3.c sqlite3.h
	CC="${CC}" CXX="${CXX}" ${CC} -std=gnu99 -c stitch.c sqlite3.c -DSTITCH_PARALLEL 
	ar rc libstitch.a stitch.o sqlite3.o

stitch_test: stitch_test.c stitch.lib
	${CC} -o stitch_test stitch_test.c -DSTITCH_PARALLEL -L. -lstitch -ldl -lpthread

clean:
	rm -f *.{o,so} libstitch.a stitch_test *.st

