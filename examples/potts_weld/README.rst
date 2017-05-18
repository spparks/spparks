==================
Running Example(s) 
==================

There are 2d and 3d examples elliptical and teardrop shaped pools.
Details of running the examples are described below for one case (3d elliptical
pool) with and without microstructure initialization.  In the examples 
the spparks was compiled into the executable *spk_flamer.gnu*.


.. _`set`: http://spparks.sandia.gov/doc/set.html
.. _`read_sites`: http://spparks.sandia.gov/doc/read_sites.html

Without base metal microstructure initialization
++++++++++++++++++++++++++++++++++++++++++++++++

To simulate grain growth during welding, without 
an initial microstructure representing base metal 
of welded material, run *spparks*:

.. code-block:: bash

   mpiexec -np 16 spk_flamer.gnu < in.steady_elliptical_weld

Use a tool to view the images: top.*.jpg


With an initial equiaxed microstructure
+++++++++++++++++++++++++++++++++++++++

In this case, run the Potts model first:

.. code-block:: bash

   mpiexec -np 16 spk_flamer.gnu < in.potts_init


Using the example input file *in.potts_init*, *spparks* will dump the file
*potts_init.dump* in the above run.  Process the dump file to write another file 
called *site.init*.  Run the script (no args).

.. code-block:: bash

   write_site_init.sh potts_init.dump site.init


Edit *in.steady_elliptical_weld* and comment out the 
`set`_ command and uncomment out the `read_sites`_ 
command; then run the potts/weld app:

.. code-block:: bash

   mpiexec -np 16 spk_flamer.gnu < in.steady_elliptical_weld

Use a tool to view the images: top.*.jpg; note the 
differences in microstructures -- especially in the 
heat affected zone.
