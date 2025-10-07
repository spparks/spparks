import sqlite3
import numpy
import struct

import stitch.libstitch as libstitch

conn = sqlite3.connect ('test.st')
c = conn.cursor ()
for row in c.execute ('select x_min, x_max, y_min, y_max, z_min, z_max, timestamp, state from blocks order by timestamp, x_min, y_min, z_min'):
    state = [x for x in struct.unpack ('iiiiiiiiiiii',row[7])]
    print ('(',row[0],',',row[2],',',row[4],')-(',row[1],',',row[3],',',row[5],'):',row[6], ' ', state)

c.close ()

(rc, f) = libstitch.open ('test.st')
(rc, field_id) = libstitch.query_field (f, 'spin')
b1=numpy.fromiter([0,2,0,2,0,3],dtype=numpy.int32).reshape(2,3,order='F')
(rc, s1, new_time) = libstitch.read_block (f, field_id, 3.0, b1)
print (b1)
print (s1)
libstitch.close (f)
