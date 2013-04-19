#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generates the user-friendly band-merged 'iphasSource' catalogues.

This script will merge the same-epoch H-alpha/r/i exposures of each IPHAS
field into band-merged catalogues, following the UKIDSS format.

This script also applies the global calibration to the magnitudes!

TODO:
 - add htmID?
"""
from __future__ import print_function, division
import os
import sys
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

HOSTNAME = os.uname()[1]
if HOSTNAME == 'uhppc11.herts.ac.uk':  # testing
    # Where are the pipeline-reduced catalogues?
    DATADIR = "/home/gb/tmp/iphas-dr2/detected"
    # Where to write the output catalogues?
    DESTINATION = "/home/gb/tmp/iphas-dr2/bandmerged"
else:  # production
    DATADIR = "/car-data/gb/iphas-dr2/detected"
    DESTINATION = "/car-data/gb/iphas-dr2/bandmerged"

# Dir of this script
SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))

# How to execute stilts?
STILTS = 'nice java -Xmx500M -XX:+UseConcMarkSweepGC -jar {0}'.format(
                                 os.path.join(SCRIPTDIR, '../lib/stilts.jar'))

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

        config = {'stilts': STILTS,
                  'runr': self.get_catalogue_path(self.run_r),
                  'runi': self.get_catalogue_path(self.run_i),
                  'runha': self.get_catalogue_path(self.run_ha),
                  'fieldid': self.fieldid,
                  'shift_r': self.shift_r,
                  'shift_i': self.shift_i,
                  'shift_ha': self.shift_ha,
                  'ocmd': os.path.join(SCRIPTDIR, 'stilts-band-merging.cmd'),
                  'output': self.output}

        cmd = """{stilts} tmatchn matcher=sky params=0.5 multimode=group \
                  nin=3 in1={runr} in2={runi} in3={runha} \
                  join1=always join2=always join3=always \
                  values1='ra dec' values2='ra dec' values3='ra dec' \
                  icmd1='setparam fieldID "{fieldid}";
                         select "aperMag2Err > 0 & aperMag2Err < 1
                                 & aperMag3Err > 0 & aperMag3Err < 1";
                         replacecol peakMag  "toFloat(peakMag  + {shift_r})";
                         replacecol aperMag2 "toFloat(aperMag2 + {shift_r})";
                         replacecol aperMag3 "toFloat(aperMag3 + {shift_r})";' \
                  icmd2='select "aperMag2Err > 0 & aperMag2Err < 1
                                 & aperMag3Err > 0 & aperMag3Err < 1"; \
                         replacecol peakMag  "toFloat(peakMag  + {shift_i})";
                         replacecol aperMag2 "toFloat(aperMag2 + {shift_i})";
                         replacecol aperMag3 "toFloat(aperMag3 + {shift_i})";' \
                  icmd3='select "aperMag2Err > 0 & aperMag2Err < 1
                                 & aperMag3Err > 0 & aperMag3Err < 1"; \
                         replacecol peakMag  "toFloat(peakMag  + {shift_ha})";
                         replacecol aperMag2 "toFloat(aperMag2 + {shift_ha})";
                         replacecol aperMag3 "toFloat(aperMag3 + {shift_ha})";' \
                  ocmd=@{ocmd} \
                  progress=none \
                  out='{output}'""".format(**config)

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


def run_all(lon1=20, lon2=220, ncores=4):
    """ Band-merge all fields """
    log.info('Band-merging in longitude range [{0},{1}['.format(lon1, lon2))

    # Make sure the output directory exists
    if not os.path.exists(DESTINATION):
        os.makedirs(DESTINATION)

    # Which fields do we want to merge?
    idx = np.where(IPHASQC.field('is_pdr')
                   & (IPHASQC.field('l') >= lon1)
                   & (IPHASQC.field('l') < lon2))[0]
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

    # Which longitude strip to process? (optional)
    if len(sys.argv) > 1:
        lon1 = int(sys.argv[1])
        lon2 = lon1 + 10
        # There are a tiny bit of sources near the survey edge at l ~ 29
        # so we add them to the strip at 30
        if lon1 == 30:
            lon1 = 25.0
    else:
        lon1 = 0
        lon2 = 360

    if HOSTNAME == 'uhppc11.herts.ac.uk':  # testing
        #run_one('5089o_jun2005')
        #run_one('3561_nov2003')
        run_all(lon1=lon1, lon2=lon2, ncores=4)
        #run_all(lon1=208, lon2=209, ncores=6)
    else:  # production
        run_all(lon1=lon1, lon2=lon2, ncores=8)
