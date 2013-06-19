#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Adds H-alpha overrides to existing zeropoint overrides

i.e. enforces zp(r) - zp(Halpha) = 3.14
"""
from __future__ import division, print_function, unicode_literals
from astropy.io import ascii
import numpy as np
import os
import sys

sys.path.append('/home/gb/dev/iphas-dr2')
from dr2 import constants

filename_target = os.path.join(constants.PACKAGEDIR, 'lib',
                               'zeropoint-overrides.csv')

filename_runs = os.path.join(constants.DESTINATION, 'runs.csv')
runs = ascii.read(filename_runs)

out = file(filename_target, 'w')
out.write('run,zp\n')

# Override each H-alpha zeropoint by enforcing zp(r) - zp(Halpha) = 3.14
for row in runs:
    if row['WFFBAND'] == 'Halpha':
        # Figure out the index of r-band runs in the same night
        idx_r = np.argwhere((np.abs(row['MJD-OBS'] - runs.field('MJD-OBS')) < 0.4)
                            & (runs.field('WFFBAND') == 'r'))[0][0]

        zp = runs[idx_r]['MAGZPT'] - 3.14
        out.write('{0},{1}\n'.format(row['run'], zp))


# Manual overrides
zp = {}
zp[381709] = 24.64 - 0.047  # 2925o_dec2003 r
zp[381710] = 23.96 - 0.027  # 2925o_dec2003 i
zp[597862] = 24.58 + 0.041  # 2832_dec2007 r
zp[597863] = 23.93 + 0.011  # 2832_dec2007 i
zp[381679] = 24.64 - 0.039  # 2914o_dec2003 r
zp[381680] = 23.96 + 0.001  # 2914o_dec2003 i
zp[598691] = 24.59 + 0.070  # 2327_dec2007 r
zp[598692] = 23.93 + 0.015  # 2327_dec2007 i
zp[528580] = 24.33 - 0.240  # 2426_oct2006 r
zp[528581] = 23.74 - 0.411  # 2426_oct2006 i
zp[471736] = 24.24 + 0.273  # 6745_sep2005 r
zp[471737] = 23.85 - 0.091  # 6745_sep2005 i
zp[948377] = 24.69 + 0.038  # 2798_nov2012 r
zp[948378] = 24.06 + 0.050  # 2798_nov2012 i
zp[530707] = 24.49 + 0.044  # 3367o_oct2006 r
zp[530708] = 23.85 + 0.038  # 3367o_oct2006 i
zp[530659] = 24.49 + 0.049  # 3352o_oct2006 r
zp[530660] = 23.85 + 0.032  # 3352o_oct2006 i
zp[430347] = 24.41 - 0.017  # 1385o_oct2004 r
zp[430348] = 23.71 + 0.060  # 1385o_oct2004 i
zp[486267] = 24.49 + 0.074  # 2694o_dec2005 r
zp[486268] = 23.84 - 0.003  # 2694o_dec2005 i
zp[381202] = 24.64 - 0.064  # 2817o_dec2003 r
zp[381203] = 23.96 - 0.036  # 2817o_dec2003 i
zp[486267] = 24.49 + 0.074  # 2694o_dec2005 r
zp[486268] = 23.84 - 0.003  # 2694o_dec2005 i
zp[381256] = 24.64 - 0.106  # 2831o_dec2003 r
zp[381257] = 23.96 - 0.097  # 2831o_dec2003 i

for myrun in zp.keys():
    out.write('{0},{1}\n'.format(myrun, zp[myrun]))

out.close()
