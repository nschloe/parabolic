# -*- coding: utf-8 -*-
#
from time_steppers import *

from parabolic.__about__ import (
    __version__,
    __author__,
    __author_email__
    )

import pipdated
if pipdated.needs_checking('parabolic'):
    msg = pipdated.check('parabolic', __version__)
    if msg:
        print(msg)
