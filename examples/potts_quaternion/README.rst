==================
Running Example(s) 
==================

Two example input scripts are provided which highlight the use of the app
style "potts/quaternion."  This app takes the basic Potts model app and adds
quaternions to each site to emulate the effect of crystalline orientation
during grain growth in polycrystals.  The natural symmetries of different
crystal systems are used to assign disorientation angles to grain boundaries,
see the documentation of this app for more details. A Read-Shockley model 
is used to calculate low-angle grain boundary energies.

Running the provided scripts will requires SPPARKS that was built with the
STITCH library. Details of running the examples are described below. In these
examples SPPARKS was compiled into the executable *spk_mpi*.


Running example scripts
+++++++++++++++++++++++++++++++++++++++

.. code-block:: bash

   mpirun -np 4 spk_mpi < in.potts_cubic
   python plot_stitch_cut.py potts_cubic --field=site 0 100 0 100

   mpirun -np 4 spk_mpi < in.potts_hcp_cutoff25
   python plot_stitch_cut.py potts_hcp_cutoff25 --field=site 0 100 0 100


This will produce log files and STITCH database files with labels "potts_cubic"
and "potts_hcp_cutoff25," which can be compared to similar simulations in the
examples/potts folder to see how Read-Shockley affects grain evolution.
The python commands will then visualize the content of the STITCH database files.
You will need to have the STITCH Python module installed, see the 
manual for more details.


