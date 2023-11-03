import inspect
import os
import subprocess
import sys
from functools import partial
from operator import attrgetter

REVISION_CMD = "git rev-parse --short HEAD"


def _get_git_revision():
    try:
        revision = subprocess.check_output(REVISION_CMD.split()).strip()
    except (subprocess.CalledProcessError, OSError):
        print("Failed to execute git to get revision")
        return None
    return revision.decode("utf-8")


def _linkcode_resolve(domain, info, package, url_fmt, revision):
    # Continue as normal for non-attributes
    if domain not in ("py", "pyx") or not info.get("module") or '.' not in info["fullname"]:
        return None
    
    class_name, _, attr = info["fullname"].partition('.')
    module = __import__(info["module"], fromlist=[class_name])
    
    # If we're trying to get an attribute of a class, skip it
    if attr:
        return None
    
    obj = getattr(module, class_name, None)
    if obj is None:
        return None
    
    obj = inspect.unwrap(obj)

    try:
        fn = inspect.getsourcefile(obj)
    except Exception:
        fn = None
    if not fn:
        try:
            fn = inspect.getsourcefile(sys.modules[obj.__module__])
        except Exception:
            fn = None
    if not fn:
        return

    # Don't include filenames from outside this package's tree
    if os.path.dirname(__import__(package).__file__) not in fn:
        return

    fn = os.path.relpath(
        fn, start=os.path.dirname(__import__(package).__file__)
    )
    try:
        lineno = inspect.getsourcelines(obj)[1]
    except Exception:
        lineno = ""
    return url_fmt.format(
        revision=revision, package=package, path=fn, lineno=lineno
    )


def make_linkcode_resolve(package, url_fmt):
    """Return a linkcode_resolve function for the given URL format.

    revision is a git commit reference (hash or name)

    package is the name of the root module of the package

    url_fmt is along the lines of ('https://github.com/USER/PROJECT/'
                                   'blob/{revision}/{package}/'
                                   '{path}#L{lineno}')
    """
    revision = _get_git_revision()
    return partial(
        _linkcode_resolve, revision=revision, package=package, url_fmt=url_fmt
    )
