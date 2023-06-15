try:
    from ._version import __version__  # noqa: F401
except ImportError:
    __version__ = "0+unknown"

# list all submodules available in omicslearn and version
__all__ = [
    "__version__",
]
