Build and install python module
===============================

Uses *build* module to create wheel.  Does not install dependencies (-n).
Builds and copies wheel to the "dist" directory.

.. code-block:: bash

  python -m build --wheel -n
  pip install dist/stitch-<VERSION>-<ARCH>.whl

For development, it may useful to build an editable project (-e).

.. code-block:: bash

  python -m build --wheel -n
  pip install --no-deps -e .

Build module api documentation
==============================

.. code-block:: bash

  sphinx-build -M html docs/source docs/


Stitch Python Module API
========================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   autoapi/index



.. JAK commented out this junk
.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


