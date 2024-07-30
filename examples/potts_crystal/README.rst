==================
Running Example(s) 
==================

* MEG QUESTION: this is in 'rst' format (sphinx etc.), this will remain the
same yes?

Two example input scripts are provided which highlight the use of the app style "potts/grain_growth" which can be used to add grain orientation to a SPPARKS Potts simulation. This app takes the basic Potts model app and adds quaternions to each grain to support emulation of grain growth. The natural symmetries of different crystal systems are used to assign disorientation values to grain boundaries. 
There are currently two variants implemented: one representing cubic systems (face-centered cubic, body-centered cubic, and simple cubic), and one for hexagonal close-packed systems. 
Running the provided scripts will require the STITCH library. Details of running the examples are described below. 
In these examples SPPARKS was compiled into the executable *spk_mpi*.

*** MEG START HERE NEXT TIME ***

Running in.potts_*_orientation scripts
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
