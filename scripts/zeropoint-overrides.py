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

filename_brent = os.path.join(constants.PACKAGEDIR, 'lib',
                              'zeropoint-overrides-brent.csv')
filename_target = os.path.join(constants.PACKAGEDIR, 'lib',
                               'zeropoint-overrides.csv')
existing_overrides = ascii.read(filename_brent)


filename_runs = os.path.join(constants.DESTINATION, 'runs.csv')
filename_runs = '/home/gb/tmp/test.csv'
runs = ascii.read(filename_runs)

os.system('cp {0} {1}'.format(filename_brent, filename_target))
out = file(filename_target, 'a')

# Override each H-alpha zeropoint by enforcing zp(r) - zp(Halpha) = 3.14
for row in runs:
    if row['WFFBAND'] == 'Halpha':
        # Ignore existing overrides
        if row['run'] in existing_overrides.field('run'):
            pass

        # Figure out the index of r-band runs in the same night
        idx_r = np.argwhere((np.abs(row['MJD-OBS'] - runs.field('MJD-OBS')) < 0.2)
                            & (runs.field('WFFBAND') == 'r'))[0][0]

        zp = runs[idx_r]['MAGZPT'] - 3.14
        out.write('{0},{1}\n'.format(row['run'], zp))

out.close()
