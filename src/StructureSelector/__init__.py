
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


from . import diff_map
from . import prepper


__all__ = ['diff_map', 'prepper']
