# Make.sh = update  style_*.h files
# Syntax: sh Make.sh style

# function to create one style_*.h file
# must whack *.d files that depend on style_*.h file,
# else Make will not recreate them

style () {
  list=`grep -l $1 $2*.h`
  if (test -e style_$3.tmp) then
    rm -f style_$3.tmp
  fi
  for file in $list; do
    qfile="\"$file\""
    echo "#include $qfile" >> style_$3.tmp
  done
  if (test ! -e style_$3.tmp) then
    rm -f style_$3.h
    touch style_$3.h
  elif (test ! -e style_$3.h) then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
    rm -f Obj_*/lammps.d
  elif (test "`diff --brief style_$3.h style_$3.tmp`" != "") then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
    rm -f Obj_*/lammps.d
  else
    rm -f style_$3.tmp
  fi
}

# create individual style files

if (test $1 = "style") then

  style APP_CLASS       app_        app        input
  style COMMAND_CLASS   ""          command    input
  style DIAG_CLASS      diag_       diag       input
  style DUMP_CLASS      dump_       dump       output
  style PAIR_CLASS      pair_       pair       potential
  style REGION_CLASS    region_     region     domain
  style SOLVE_CLASS     solve_      solve      input

fi
