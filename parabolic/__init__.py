# -*- coding: utf-8 -*-
#
from __future__ import print_function

# pylint: disable=wildcard-import
from .time_steppers import *

from .__about__ import (
    __version__,
    __author__,
    __author_email__
    )

try:
    import pipdate
except ImportError:
    pass
else:
    if pipdate.needs_checking(__name__):
        print(pipdate.check(__name__, __version__), end='')
