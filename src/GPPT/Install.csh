# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_gppt.h ..

  cp app_gppt.cpp ..
  cp fitness.cpp ..
  cp population.cpp ..
  cp solve_gppt.cpp ..

  cp app_gppt.h ..
  cp fitness.h ..
  cp population.h ..
  cp solve_gppt.h ..

else if ($1 == 0) then

  rm ../style_gppt.h
  touch ../style_gppt.h

  rm ../app_gppt.cpp
  rm ../fitness.cpp
  rm ../population.cpp
  rm ../solve_gppt.cpp

  rm ../app_gppt.h
  rm ../fitness.h
  rm ../population.h
  rm ../solve_gppt.h

endif
