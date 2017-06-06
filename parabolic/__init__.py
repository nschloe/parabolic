# -*- coding: utf-8 -*-
#
from __future__ import print_function

import pipdated

# pylint: disable=wildcard-import
from .time_steppers import *

from .__about__ import (
    __version__,
    __author__,
    __author_email__
    )

if pipdated.needs_checking(__name__):
    print(pipdated.check(__name__, __version__))
