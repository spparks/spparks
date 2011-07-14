# Make.csh = update Makefile.lib or Makefile.list or style_*.h files
# use current list of *.cpp and *.h files in src dir

if ($1 == "Makefile.lib") then

  set list1 = `ls [abcde]*.cpp`
  set list2 = `ls [fghij]*.cpp`
  set list3 = `ls [klmno]*.cpp | sed s/^main\.cpp//`
  set list4 = `ls [pqrst]*.cpp`
  set list5 = `ls [uvwxyz]*.cpp`
  sed -i -e "s/SRC =	.*/SRC =	$list1/" Makefile.lib
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list2/" Makefile.lib
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list3/" Makefile.lib
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list4/" Makefile.lib
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list5/" Makefile.lib

  set list1 = `ls [abcde]*.h`
  set list2 = `ls [fghij]*.h`
  set list3 = `ls [klmno]*.h`
  set list4 = `ls [pqrst]*.h`
  set list5 = `ls [uvwxyz]*.h`
  sed -i -e "s/INC =	.*/INC =	$list1/" Makefile.lib
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list2/" Makefile.lib
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list3/" Makefile.lib
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list4/" Makefile.lib
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list5/" Makefile.lib

else if ($1 == "Makefile.list") then

  set list1 = `ls [abcde]*.cpp`
  set list2 = `ls [fghij]*.cpp`
  set list3 = `ls [klmno]*.cpp`
  set list4 = `ls [pqrst]*.cpp`
  set list5 = `ls [uvwxyz]*.cpp`
  sed -i -e "s/SRC =	.*/SRC =	$list1/" Makefile.list
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list2/" Makefile.list
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list3/" Makefile.list
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list4/" Makefile.list
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list5/" Makefile.list

  set list1 = `ls [abcde]*.h`
  set list2 = `ls [fghij]*.h`
  set list3 = `ls [klmno]*.h`
  set list4 = `ls [pqrst]*.h`
  set list5 = `ls [uvwxyz]*.h`
  sed -i -e "s/INC =	.*/INC =	$list1/" Makefile.list
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list2/" Makefile.list
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list3/" Makefile.list
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list4/" Makefile.list
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list5/" Makefile.list

else if ($1 == "style") then

  set list = `grep -l APP_CLASS app_*.h`
  if (-e style_app.tmp) then
    rm style_app.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_app.tmp
  end
  if (! -e style_app.tmp) then
     rm style_app.h
     touch style_app.h
  else if (! -e style_app.h) then
     mv style_app.tmp style_app.h
     rm Obj_*/input.d
  else if (`diff style_app.h style_app.tmp` != "") then
     mv style_app.tmp style_app.h
     rm Obj_*/input.d
  else
     rm style_app.tmp
  endif

  set list = `grep -l COMMAND_CLASS *.h`
  if (-e style_command.tmp) then
    rm style_command.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_command.tmp
  end
  if (! -e style_command.tmp) then
     rm style_command.h
     touch style_command.h
  else if (! -e style_command.h) then
     mv style_command.tmp style_command.h
     rm Obj_*/input.d
  else if (`diff style_command.h style_command.tmp` != "") then
     mv style_command.tmp style_command.h
     rm Obj_*/input.d
  else
     rm style_command.tmp
  endif

  set list = `grep -l DIAG_CLASS diag_*.h`
  if (-e style_diag.tmp) then
    rm style_diag.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_diag.tmp
  end
  if (! -e style_diag.tmp) then
     rm style_diag.h
     touch style_diag.h
  else if (! -e style_diag.h) then
     mv style_diag.tmp style_diag.h
     rm Obj_*/input.d
  else if (`diff style_diag.h style_diag.tmp` != "") then
     mv style_diag.tmp style_diag.h
     rm Obj_*/input.d
  else
     rm style_diag.tmp
  endif

  set list = `grep -l DUMP_CLASS dump_*.h`
  if (-e style_dump.tmp) then
    rm style_dump.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_dump.tmp
  end
  if (! -e style_dump.h) then
     mv style_dump.tmp style_dump.h
     rm Obj_*/input.d
  else if (`diff style_dump.h style_dump.tmp` != "") then
     mv style_dump.tmp style_dump.h
     rm Obj_*/input.d
  else
     rm style_dump.tmp
  endif

  set list = `grep -l PAIR_CLASS pair_*.h`
  if (-e style_pair.tmp) then
    rm style_pair.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_pair.tmp
  end
  if (! -e style_pair.tmp) then
     rm style_pair.h
     touch style_pair.h
  else if (! -e style_pair.h) then
     mv style_pair.tmp style_pair.h
     rm Obj_*/potential.d
  else if (`diff style_pair.h style_pair.tmp` != "") then
     mv style_pair.tmp style_pair.h
     rm Obj_*/potential.d
  else
     rm style_pair.tmp
  endif

  set list = `grep -l REGION_CLASS region_*.h`
  if (-e style_region.tmp) then
    rm style_region.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_region.tmp
  end
  if (! -e style_region.tmp) then
     rm style_region.h
     touch style_region.h
  else if (! -e style_region.h) then
     mv style_region.tmp style_region.h
     rm Obj_*/domain.d
  else if (`diff style_region.h style_region.tmp` != "") then
     mv style_region.tmp style_region.h
     rm Obj_*/domain.d
  else
     rm style_region.tmp
  endif

  set list = `grep -l SOLVE_CLASS solve_*.h`
  if (-e style_solve.tmp) then
    rm style_solve.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_solve.tmp
  end
  if (! -e style_solve.tmp) then
     rm style_solve.h
     touch style_solve.h
  else if (! -e style_solve.h) then
     mv style_solve.tmp style_solve.h
     rm Obj_*/input.d
  else if (`diff style_solve.h style_solve.tmp` != "") then
     mv style_solve.tmp style_solve.h
     rm Obj_*/input.d
  else
     rm style_solve.tmp
  endif

endif
