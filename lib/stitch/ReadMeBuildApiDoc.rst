= Howto Build Stitch Api Documentation using SPHINX-Build

Command below builds the documentation.  Location of *sphinx* related files describing 
stitch is docs/source. *stitch* api is created into the directory docs/html.  

libstitch api documentation is in src/stitch/libstitch/libstitch.pyi

While it may be possible to include the documentation with the source code directly, 
this approach was not used because tools are not quite there for documenting 
a native c-api.  Instead the python api documentation is written in 
the file stitch/libstitch/libstitch.pyi.  The default build of the libstitch api 
documentation looks standard with the *sphinx* default look and feel.

Using sphinx-build to build api documentation.  

. code-block:: bash

  sphinx-build -M html docs/source docs/html

Important files for the python stitch api documentation.

* pyproject.toml
* src/stitch/libstitch.pyi
* docs/source/conf.py
* docs/source/index.rst

== TODO Identify dependencies for build api documentation
Identify all dependencies for building the python api documentation with
sphinx-build.  Nominally, docs/source/conf.py lists a set of dependencies but
some work is needed when installing dependencies with pip as there is not a
one-to-one correpsondence with the dependency names in conf.py and what pip
installs.  Its not hard to figure out.

