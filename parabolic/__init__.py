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

if pipdated.needs_checking('parabolic'):
    msg = pipdated.check('parabolic', __version__)
    if msg:
        print(msg)
