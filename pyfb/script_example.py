#!/usr/bin/env python

import sys
import os
sys.path.append(os.environ["FIREBALLHOME"])

from pyfb.geometry.dinamic import *

din=dinamic()
din.load_xyz("answer.xyz")
din.print()


