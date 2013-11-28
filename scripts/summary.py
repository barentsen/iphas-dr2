#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Prints a summary of the contents of the IPHAS source catalogue.
"""
import os
from astropy.io import fits
from astropy import log
import numpy as np
import sys
from dr2 import constants

n_sources = 0
n_r20 = 0
n_reliable = 0
n_deblend = 0
n_reliable_deblend = 0
n_pair = 0
n_saturated = 0
n_brightNeighb = 0

path = os.path.join(constants.DESTINATION, 'concatenated', 'full')
for filename in os.listdir(path):
    if filename.endswith('fits.gz'):
        print filename
        myfile = os.path.join(path, filename)
        log.info(myfile)
        f = fits.open(myfile)
        #n_sources += f[1].header['NAXIS2']
        n_sources += f[1].data['ra'].size
        n_r20 += (f[1].data['r'] < 21).sum()
        n_reliable += f[1].data['reliable'].sum()
        n_reliable_deblend += (f[1].data['reliable'] & f[1].data['deblend']).sum()
        n_deblend += f[1].data['deblend'].sum()
        n_pair += (f[1].data['sourceID2'] != ' ').sum()
        n_saturated += f[1].data['saturated'].sum()
        n_brightNeighb += f[1].data['brightNeighb'].sum()
        print "{0} sources so far".format(n_sources)

with open('summary.txt', 'w') as out:
    out.write("#Unique sources: {0}\n".format(n_sources))
    out.write("#Sources r < 21: {0}\n".format(n_r20))
    out.write("#Reliable sources: {0}\n".format(n_reliable))
    out.write("#Deblend sources: {0}\n".format(n_deblend))
    out.write("#Reliable+deblend: {0}\n".format(n_reliable_deblend))
    out.write("#Paired sources: {0}\n".format(n_pair))
    out.write("#Saturated sources: {0}\n".format(n_saturated))
    out.write("#Bright neighb sources: {0}\n".format(n_brightNeighb))
 
