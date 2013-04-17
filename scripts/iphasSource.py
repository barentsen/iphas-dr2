#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generates the user-friendly band-merged 'iphasSource' catalogues.

This script will merge the same-epoch H-alpha/r/i exposures of each IPHAS
field into band-merged catalogues, following the UKIDSS format.

This script also applies the global calibration to the magnitudes!

TODO:
 - add priOrSec;
 - add htmID?
"""
from __future__ import print_function, division
import os
import numpy as np
from astropy.io import fits
from astropy import log
from multiprocessing import Pool

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Hywel Farnhill', 'Robert Greimel', 'Janet Drew',
               'Cambridge Astronomical Survey Unit']


#############################
# CONSTANTS & CONFIGURATION
#############################

if os.uname()[1] == 'uhppc11.herts.ac.uk':
    # Where are the pipeline-reduced catalogues?
    DATADIR = "/home/gb/tmp/iphas-dr2/iphasDetection"
    # Where to write the output catalogues?
    DESTINATION = "/home/gb/tmp/iphas-dr2/iphasSource"
else:
    DATADIR = "/car-data/gb/iphas-dr2/iphasDetection"
    DESTINATION = "/car-data/gb/iphas-dr2/iphasSource"

# How to execute stilts?
STILTS = 'nice java -Xmx500M -XX:+UseConcMarkSweepGC -jar lib/stilts.jar'

# Where is the IPHAS quality control table?
IPHASQC = fits.getdata('/home/gb/dev/iphas-qc/qcdata/iphas-qc.fits', 1)


###########
# CLASSES
###########

class IPHASException(Exception):
    """Raised when a field cannot be band-merged."""
    pass


class BandMerge():
    """
    Class to read in iphasSource tables and perform a band-merge.
    """

    def __init__(self, fieldid, run_ha, run_r, run_i,
                 shift_ha=0.0, shift_r=0.0, shift_i=0.0):
        """Constructor"""
        self.fieldid = fieldid
        self.run_ha = run_ha
        self.run_r = run_r
        self.run_i = run_i
        self.shift_ha = shift_ha
        self.shift_r = shift_r
        self.shift_i = shift_i
        self.output = os.path.join(DESTINATION, fieldid+'.fits')

    def get_catalogue_path(self, run):
        """Returns the full path of a run's 'iphasDetection' catalogue."""
        return os.path.join(DATADIR, '{}_det.fits'.format(run))

    def get_stilts_command(self):
        """Returns the stilts command used to perform a band-merge."""
        cmd = """{0} tmatchn matcher=sky params=0.5 multimode=group nin=3 \
                  in1={1} in2={2} in3={3} \
                  join1=always join2=always join3=always \
                  values1='ra dec' values2='ra dec' values3='ra dec' \
                  icmd1='setparam fieldID "{4}";
                         select "aperMag2Err > 0 & aperMag2Err < 1
                                 & aperMag3Err > 0 & aperMag3Err < 1";
                         replacecol peakMag  "toFloat(peakMag  + {5})";
                         replacecol aperMag2 "toFloat(aperMag2 + {5})";
                         replacecol aperMag3 "toFloat(aperMag3 + {5})";' \
                  icmd2='select "aperMag2Err > 0 & aperMag2Err < 1
                                 & aperMag3Err > 0 & aperMag3Err < 1"; \
                         replacecol peakMag  "toFloat(peakMag  + {6})";
                         replacecol aperMag2 "toFloat(aperMag2 + {6})";
                         replacecol aperMag3 "toFloat(aperMag3 + {6})";' \
                  icmd3='select "aperMag2Err > 0 & aperMag2Err < 1
                                 & aperMag3Err > 0 & aperMag3Err < 1"; \
                         replacecol peakMag  "toFloat(peakMag  + {7})";
                         replacecol aperMag2 "toFloat(aperMag2 + {7})";
                         replacecol aperMag3 "toFloat(aperMag3 + {7})";' \
                  ocmd=@stilts-band-merging.cmd \
                  progress=none \
                  out='{8}'""".format(STILTS,
                                      self.get_catalogue_path(self.run_r),
                                      self.get_catalogue_path(self.run_i),
                                      self.get_catalogue_path(self.run_ha),
                                      self.fieldid,
                                      self.shift_r,
                                      self.shift_i,
                                      self.shift_ha,
                                      self.output)
        return cmd

    def run(self):
        """Perform the band-merging."""
        cmd = self.get_stilts_command()
        status = os.system(cmd)
        return status


###########
# FUNCTIONS
###########


def run_one(fieldid):
    """ Band-merge a single field """
    # Which index does the field have in the QC table?
    idx = np.where(IPHASQC.field('id') == fieldid)[0]
    if len(idx) != 1:
        raise IPHASException('{}: error identifying runs'.format(fieldid))

    # Which calibration shifts to apply?
    shifts = {}
    for band in ['h', 'r', 'i']:
        shifts[band] = (IPHASQC.field('zp{}_pdr'.format(band))[idx[0]]
                        - IPHASQC.field('zp{}'.format(band))[idx[0]])
        if np.isnan(shifts[band]):
            shifts[band] = 0.0

    # Carry out the band-merging
    bm = BandMerge(fieldid,
                   IPHASQC.field('run_ha')[idx[0]],
                   IPHASQC.field('run_r')[idx[0]],
                   IPHASQC.field('run_i')[idx[0]],
                   shifts['h'], shifts['r'], shifts['i'])
    status = bm.run()

    log.info('{}: {}'.format(fieldid, status))
    return status


def run_all(ncores=4):
    """ Band-merge all fields """
    # Which fields do we want to merge?
    idx = np.where(IPHASQC.field('is_pdr'))[0]
    fieldids = IPHASQC.field('id')[idx]
    # Distribute the work over ncores
    p = Pool(processes=ncores)
    results = p.imap(run_one, fieldids)
    for i in results:
        pass


###################
# MAIN EXECUTION
###################

if __name__ == '__main__':
    run_one('5089o_jun2005')
    #run_all(8)
