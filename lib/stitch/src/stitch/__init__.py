from importlib import metadata

try:
    __version__=metadata.version(__package__)

except metadata.PackageNotFoundError:
    __version__ = "0.0.0+unknown"

del metadata
