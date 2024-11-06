==================
Running Example(s) 
==================

Four example input scripts are provided which highlight the use of the app style "potts/quaternion."  This app takes the basic Potts model app and adds quaternions to each site to emulate the effect of crystalline orientation during grain growth in polydrystals.  The natural symmetries of different crystal systems are used to assign disorientation angles to grain boundaries, see the documentation of this app for more details.


There are currently two types of file included for each crystal symmetry type: one file that uses the Read-Shockley model to calculate low-angle grain boundary energies, and one file that does not for comparison.  


Running the provided scripts will require the STITCH library. Details of running the examples are described below. 
In these examples SPPARKS was compiled into the executable *spk_mpi*.


Running example scripts
+++++++++++++++++++++++++++++++++++++++

.. code-block:: bash

   mpirun -np 4 spk_mpi < in.potts_cubic
   mpirun -np 4 spk_mpi < in.potts_hcp_cutoff25


This will produce log files and STITCH database files with labels "potts_cubic" and "potts_hcp_cutoff25," which can be compared to similar simulations in the examples/potts folder to see how Read-Shockley affects the simulation output. 
The content of the STITCH database files can be visualized using the STITCH Python module, see e.g. examples/potts_am_bezier/plot_stitch_cut.py.
