# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_gppt.h ..

  cp app_gppt.cpp ..
  cp const_node.cpp ..
  cp divide_node.cpp ..
  cp fitness.cpp ..
  cp minus_node.cpp ..
  cp node.cpp ..
  cp plus_node.cpp ..
  cp population.cpp ..
  cp pow_node.cpp ..
  cp solve_gppt.cpp ..
  cp times_node.cpp ..
  cp tree.cpp ..
  cp var_node.cpp ..

  cp app_gppt.h ..
  cp const_node.h ..
  cp divide_node.h ..
  cp fitness.h ..
  cp minus_node.h ..
  cp node.h ..
  cp plus_node.h ..
  cp population.h ..
  cp pow_node.h ..
  cp solve_gppt.h ..
  cp times_node.h ..
  cp tree.h ..
  cp var_node.h ..

else if ($1 == 0) then

  rm ../style_gppt.h
  touch ../style_gppt.h

  rm ../app_gppt.cpp
  rm ../const_node.cpp
  rm ../divide_node.cpp
  rm ../fitness.cpp
  rm ../minus_node.cpp
  rm ../node.cpp
  rm ../plus_node.cpp
  rm ../population.cpp
  rm ../pow_node.cpp
  rm ../solve_gppt.cpp
  rm ../times_node.cpp
  rm ../tree.cpp
  rm ../var_node.cpp

  rm ../app_gppt.h
  rm ../const_node.h
  rm ../divide_node.h
  rm ../fitness.h
  rm ../minus_node.h
  rm ../node.h
  rm ../plus_node.h
  rm ../population.h
  rm ../pow_node.h
  rm ../solve_gppt.h
  rm ../times_node.h
  rm ../tree.h
  rm ../var_node.h

endif
