#!/usr/bin/env python

import sys
import datetime

for line in sys.stdin:
    print(f"{datetime.datetime.now()}, {line.rstrip()}")
