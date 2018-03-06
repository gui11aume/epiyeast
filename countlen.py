#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
from gzopen import gzopen

avg = 0.0
tot = 0

with gzopen(sys.argv[1]) as f:
   for line in f:
      items = line.split()
      avg += len(items[0])
      tot += 1
print avg/tot
