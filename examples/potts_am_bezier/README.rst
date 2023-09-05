==================
Running Example(s) 
==================

Two example input scripts are provided which highlight the use of the
app style "potts/am/bezier" which can apply an experimentally fitted
bezier melt-pool profile to drive the SPPARKS simulations. This gives
the user the unique ability to incorporate experimentally parameterized
melt pool shapes into the SPPARKS simulation. Running the provided 
scripts will require the STITCH library. Details of running  the 
examples are described below. In these examples spparks was compiled 
into the executable *spk_mpi*.

Running in.potts_am_bezier_* scripts
+++++++++++++++++++++++++++++++++++++++

.. code-block:: bash

   mpiexec -np 1 spk_mpi < in.potts_am_bezier_small2d

.. code-block:: bash

   mpiexec -np 16 spk_mpi < in.potts_am_bezier_large2d

This will produce the stitch file "potts_am_bezier_small2d.st" which
can be examined using the python script "plot_stitch_cut.py". Running
the "plot_stitch_cut.py" script, as shown below, will produce a png
rendering of the SPPARKS microstructure.

.. code-block:: bash

   python3 plot_stitch_cut.py potts_am_bezier_small2d --field=site 0 100 0 100

.. code-block:: bash

   python3 plot_stitch_cut.py potts_am_bezier_large2d --field=site 0 500 0 500

An Ovito rendering can also be produced by uncommenting the following line in
each of the input files: "dump 1 text 5.0 dump.additive4.* id i1 d1 x y z"
This will produce a series of dump files "dump.additive4.*" which can be 
directly visualized in Ovito.
