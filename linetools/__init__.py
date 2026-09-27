"""linetools: analysis of 1d astronomical spectra."""

from importlib.metadata import PackageNotFoundError, version as _version

try:
    __version__ = _version("linetools")
except PackageNotFoundError:  # not installed, e.g. run from a source tree
    __version__ = "unknown"

__all__ = ["__version__"]
