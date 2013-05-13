#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Constants used in the IPHAS Data Release modules."""
import os
from astropy.io import fits

DEBUGMODE = False

RAWDATADIR = '/car-data/gb/iphas'  # Where are the pipeline-reduced catalogues?
DESTINATION = '/car-data/gb/iphas-dr2'  # Where to write output catalogues?
"""
HOSTNAME = os.uname()[1]
if HOSTNAME == 'uhppc11.herts.ac.uk':  # testing machine
    DEBUGMODE = True
    RAWDATADIR = '/media/0133d764-0bfe-4007-a9cc-a7b1f61c4d1d/iphas'
    DESTINATION = '/home/gb/tmp/iphas-dr2'
if HOSTNAME == 'gvm':  # testing machine
    DEBUGMODE = True
    RAWDATADIR = '/media/uh/media/0133d764-0bfe-4007-a9cc-a7b1f61c4d1d/iphas'
    DESTINATION = '/home/gb/tmp/iphas-dr2'
"""
PACKAGEDIR = os.path.dirname(os.path.abspath(__file__))

# Where is the IPHAS quality control table?
IPHASQC = fits.getdata('/home/gb/dev/iphas-qc/qcdata/iphas-qc.fits', 1)
IPHASQC_COND_RELEASE = (IPHASQC['is_pdr'] & (IPHASQC['qflag'] != 'D'))

# How to execute stilts?
STILTS = 'nice java -Xmx2000M -XX:+UseConcMarkSweepGC -jar {0}'.format(
                                os.path.join(PACKAGEDIR, 'lib/stilts.jar'))

# Fields within this radius will be considered to overlap
FIELD_MAXDIST = 0.8  # degrees

# Width of the Galactic Plane strip to process
STRIPWIDTH = 5  # degrees galactic longitude

# Detections within this radius will be considered identical
MATCHING_DISTANCE = 0.5  # arcsec

 #INT/WFC CCD pixel scale
PXSCALE = 0.333  # arcsec/pix

# Filter names
BANDS = ['r', 'i', 'ha']

# Which are the possible filenames of the confidence maps?
CONF_NAMES = {'Halpha': ['Ha_conf.fits', 'Ha_conf.fit',
                         'Halpha_conf.fit',
                         'ha_conf.fits', 'ha_conf.fit',
                         'h_conf.fits', 'h_conf.fit',
                         'Halpha:197_iphas_aug2003_cpm.fit',
                         'Halpha:197_iphas_sep2003_cpm.fit',
                         'Halpha:197_iphas_oct2003_cpm.fit',
                         'Halpha:197_iphas_nov2003_cpm.fit',
                         'Halpha:197_nov2003b_cpm.fit',
                         'Halpha:197_dec2003_cpm.fit',
                         'Halpha:197_jun2004_cpm.fit',
                         'Halpha:197_iphas_jul2004a_cpm.fit',
                         'Halpha:197_iphas_jul2004_cpm.fit',
                         'Halpha:197_iphas_aug2004a_cpm.fit',
                         'Halpha:197_iphas_aug2004b_cpm.fit',
                         'Halpha:197_iphas_dec2004b_cpm.fit'],
              'r': ['r_conf.fit', 'r_conf.fits',
                    'r:214_iphas_aug2003_cpm.fit',
                    'r:214_dec2003_cpm.fit',
                    'r:214_iphas_nov2003_cpm.fit',
                    'r:214_nov2003b_cpm.fit',
                    'r:214_iphas_sep2003_cpm.fit',
                    'r:214_iphas_aug2004a_cpm.fit',
                    'r:214_iphas_aug2004b_cpm.fit',
                    'r:214_iphas_jul2004a_cpm.fit',
                    'r:214_iphas_jul2004_cpm.fit',
                    'r:214_jun2004_cpm.fit'],
              'i': ['i_conf.fit', 'i_conf.fits',
                    'i:215_iphas_aug2003_cpm.fit',
                    'i:215_dec2003_cpm.fit',
                    'i:215_iphas_nov2003_cpm.fit',
                    'i:215_nov2003b_cpm.fit',
                    'i:215_iphas_sep2003_cpm.fit',
                    'i:215_iphas_aug2004a_cpm.fit',
                    'i:215_iphas_aug2004b_cpm.fit',
                    'i:215_iphas_jul2004a_cpm.fit',
                    'i:215_iphas_jul2004_cpm.fit',
                    'i:215_jun2004_cpm.fit']}
