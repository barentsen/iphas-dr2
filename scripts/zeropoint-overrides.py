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

out.close()
