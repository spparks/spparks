#!/bin/bash
# ----------------------------------------------------------------------
# Jak's backup utility
unset PATH	# avoid accidental use of $PATH

# ------------- system commands used by this script --------------------
ID=/usr/bin/id;
ECHO=/bin/echo;
MOUNT=/bin/mount;
UMOUNT=/bin/umount;
RM=/bin/rm;
MV=/bin/mv;
CP=/bin/cp;
TOUCH=/bin/touch;
RSYNC=/usr/bin/rsync;

# ------source and destination info --------------------------------------
SPPARKS=$HOME/jaks.git/spparks.devel

# Files for building libstitch.a; also source code for python module
# Stitch Source Files
DEST=$SPPARKS/lib/stitch/libstitch
SRC=$HOME/jaks.git/stitch.cloned/libstitch
# NOTE: 'setup.py' was slightly tweaked for 'SPPARKS' to exclude 
#   weld module directory; in future it may be good to add this 
#   back in here.
FILES="sqlite3.h sqlite3.c stitch.c stitch.h stitchmodule.c "
FILES+="stitch_test.c stitch_test_data.h setup.py"
for f in $FILES; do
   # DRY-RUN '-n'
   $RSYNC -va $SRC/$f $DEST ;
   #$RSYNC -n -va $SRC/$f $DEST ;
done

# Files for python module build
DEST=$SPPARKS/lib/stitch
# Python unit test
SRC=$HOME/jaks.git/stitch.cloned/
FILES="LICENSE setup.py setup.cee.cfg setup.flamer.cfg"
for f in $FILES; do
   # DRY-RUN '-n'
   echo "rsync file: " $f
   $RSYNC -va $SRC/$f $DEST ;
   #$RSYNC -n -va $SRC/$f $DEST ;
done

# Python unit test
DEST=$SPPARKS/lib/stitch/verify
SRC=$HOME/jaks.git/stitch.cloned/integrated_test/weld/potts
FILES="unit_cv_readwrite.py"
for f in $FILES; do
   # DRY-RUN '-n'
   $RSYNC -va $SRC/$f $DEST ;
   #$RSYNC -n -va $SRC/$f $DEST ;
done
