#!/usr/bin/env python3

import sys
import datetime

for line in sys.stdin:
    print(f"{datetime.datetime.now()}, {line.rstrip()}")
