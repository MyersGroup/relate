#!/usr/bin/env python

import numpy as np
import sys

filename=sys.argv[1]
f = open(filename, "r")
a = np.fromfile(f, dtype=np.int32)
print a[0], a[1], a[2]
