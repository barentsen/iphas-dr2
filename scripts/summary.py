#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Prints a summary of the contents of the IPHAS source catalogue.
"""
import os
from astropy.io import fits
from astropy import log
import sys
sys.path.append('../')
from dr2 import constants

n_sources = 0
n_reliable = 0

path = os.path.join(constants.DESTINATION, 'catalogue', 'light')
for filename in os.listdir(path):
    if filename.endswith('fits'):
        myfile = os.path.join(path, filename)
        log.info(myfile)
        f = fits.open(myfile)
        n_sources += f[1].header['NAXIS2']
        #n_reliable += f[1].data.field('reliable').sum()
        n_reliable += f[1].data.size

print "#Unique sources: {0}".format(n_sources)
print "#Reliable sources: {0}".format(n_reliable)
