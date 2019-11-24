#!/usr/bin/env python

import sys
import os
sys.path.append(os.environ["FIREBALLHOME"])

from pyfb.format.xyz import *

load_xyz("answer.xyz")
print_xyz()


