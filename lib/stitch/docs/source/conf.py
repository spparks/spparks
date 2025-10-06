
import os
import sys

# Add the path to your C extension source
sys.path.insert(0, os.path.abspath('../../src/libstitch'))

# Sphinx extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'autoapi.extension',
    'numpydoc',
    ]

# NumPy configuration
napoleon_numpy_docstring = True

# Python is default domain
primary_domain = 'py'

html_theme = "classic"

# Configure AutoAPI
autoapi_type = 'python'
autoapi_dirs = ['../../src']
autoapi_file_patterns = ['*.pyi']

# Other Sphinx settings
project = 'Stitch Python Extension'
copyright = '2025, John Mitchell and Jay Lofstead'
author = 'John Mitchell'


