===============
Running Example 
===============

.. _`set`: http://spparks.sandia.gov/doc/set.html
.. _`read_sites`: http://spparks.sandia.gov/doc/read_sites.html

Without base metal microstructure initialization
++++++++++++++++++++++++++++++++++++++++++++++++

To simulate grain growth during welding, without 
an initial microstructure representing base metal 
of welded material (each voxel begins with a randomly 
assigned spin), run *spparks*:

.. code-block:: bash

   mpiexec -np 16 spk_mac_mpi < in.weldRandom


With an initial equiaxed microstructure
+++++++++++++++++++++++++++++++++++++++

In this case, run the Potts model first:

.. code-block:: bash

   mpiexec -np 16 spk_mac_mpi < in.potts3D


Then process the dump file to write another file 
called *site.init*.  Run the script (no args).

.. code-block:: bash

   write_site_init.sh


Then run the weld model using the in.weldRead 
input file *in.weldRead*, which uses the `read_sites`_ 
command rather than the `set`_ command to initialize
an equiaxed grain microstructure:

.. code-block:: bash

   mpiexec -np 16 spk_mac_mpi < in.weldRead
