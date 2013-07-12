#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Computes magnitude offsets between overlapping exposures.

This script carries out pair-wise crossmatching of reliable objects in
overlapping exposures. It then computes the median offset between those
magnitudes, which provides the input required to carry out a global
re-calibration.

The output is a CSV file stored at "constants.DESTINATION/offsets-{band}.csv"

The columns in the CSV file are:
run1   -- reference run number
run2   -- comparison run number
offset -- median(run1_magnitudes - run2_magnitudes)
std    -- std(run1_magnitudes - run2_magnitudes)
n      -- len(run1_magnitudes - run2_magnitudes)

Only runs which are part of the data release are included.
"""
from __future__ import division, print_function, unicode_literals
import os
import sys
from multiprocessing import Pool
import numpy as np
from astropy import log
from astropy.io import fits
import util
import constants
from constants import IPHASQC

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Hywel Farnhill', 'Janet Drew']


#############################
# CONSTANTS & CONFIGURATION
#############################

# Magnitude ranges across whcih we will compute the median offsets
MAGLIMITS = {'r': [15, 17.5], 'i': [14, 16.5], 'ha': [15, 17.5]}


###########
# CLASSES
###########

class OffsetMachine(object):
    """Computes photometric offsets between exposures."""

    def __init__(self, run):
        self.run = run
        self.data = self.get_data(run)
        self.band = self.data['band']

    def get_data(self, myrun):
        f = fits.open(self.filename(myrun))
        data = {'ra': f[1].data['ra'],
                'dec': f[1].data['dec'],
                'aperMag2': f[1].data['aperMag2'],
                'errBits': f[1].data['errBits'],
                'band': f[1].data['band'][0]}
        f.close()
        return data

    def filename(self, run):
        """Returns the full path of the detection catalogue of the run."""
        return os.path.join(constants.DESTINATION,
                            'detected',
                            '{0}_det.fits'.format(run))

    def overlap_runs(self):
        """Returns overlapping runs which are in the data release and in the same filter."""
        cond_run = (IPHASQC['run_'+str(self.band)] == self.run)
        idx = np.argwhere(cond_run)[0]
        myra = IPHASQC['ra'][idx]
        mydec = IPHASQC['dec'][idx]

        dist = util.sphere_dist(myra, mydec, IPHASQC['ra'], IPHASQC['dec'])
        idx2 = (constants.IPHASQC_COND_RELEASE
                & (dist < constants.FIELD_MAXDIST)
                & (IPHASQC['run_'+str(self.band)] != self.run))

        return IPHASQC['run_'+str(self.band)][idx2]

    def relative_offsets(self):
        """Returns the offsets."""
        log.debug('Computing offsets for run {0}'.format(self.run))
        offsets = []
        for run2 in self.overlap_runs():
            log.debug(str(run2))
            result = self._compute_relative_offsets(run2)
            if result is not None:
                offsets.append(result)
        return offsets

    def _compute_relative_offsets(self, run2):
        """Computes the offset against a specified run."""
        limit_bright = MAGLIMITS[self.band][0]
        limit_faint = MAGLIMITS[self.band][1]

        offsets = []

        offset_data = self.get_data(run2)

        cond_reliable1 = ((self.data['aperMag2'] > limit_bright)
                          & (self.data['aperMag2'] < limit_faint)
                          & (self.data['errBits'] == 0))
        cond_reliable2 = ((offset_data['aperMag2'] > limit_bright-0.2)
                          & (offset_data['aperMag2'] < limit_faint+0.2)
                          & (offset_data['errBits'] == 0))

        for idx1 in np.argwhere(cond_reliable1):
            idx2 = util.crossmatch(self.data['ra'][idx1],
                                   self.data['dec'][idx1],
                                   offset_data['ra'][cond_reliable2],
                                   offset_data['dec'][cond_reliable2])
            if idx2 is not None:
                offsets.append(self.data['aperMag2'][idx1]
                               - offset_data['aperMag2'][cond_reliable2][idx2])

        if len(offsets) < 5:
            return None
        else:
            offsets = np.array(offsets)
            return {'run1': self.run,
                    'run2': run2,
                    'offset': np.median(offsets),
                    'std': np.std(offsets),
                    'n': len(offsets)}


###########
# FUNCTIONS
###########

def offsets_relative_one(run):
    try:
        om = OffsetMachine(run)
        return om.relative_offsets()
    except Exception, e:
        log.error('UNEXPECTED EXCEPTION FOR RUN {0}: {1}'.format(run, e))
        return [None]


def offsets_relative(band, ncores=4):
    """Writes the file offsets-relative.csv with overlap offsets."""
    assert(band in constants.BANDS)
    log.info('Starting to compute offsets for band {0}'.format(band))

    # Write the results
    filename = os.path.join(constants.DESTINATION,
                            'offsets-{0}.csv'.format(band))
    out = open(filename, 'w')
    out.write('run1,run2,offset,std,n\n')

    # Distribute the work over ncores
    p = Pool(processes=ncores)
    runs = IPHASQC['run_'+str(band)][constants.IPHASQC_COND_RELEASE]
    results = p.map(offsets_relative_one, runs)
    for i, field in enumerate(results):
        log.info('Completed run {0}/{1}'.format(i, len(runs)))
        for row in field:
            if row is not None:
                out.write('{run1},{run2},{offset},{std},{n}\n'.format(**row))

        if (i % 20) == 0:  # Make sure we write at least every 20 runs
            out.flush()

    out.close()


###################
# MAIN EXECUTION
###################

if __name__ == '__main__':

    # Which band to process?
    if len(sys.argv) > 1:
        band = sys.argv[1]
    else:
        raise Exception('Please give the band as the first argument')

    #constants.DEBUGMODE = True
    if constants.DEBUGMODE:
        log.setLevel('DEBUG')
        #o = offsets_relative_one(484350)
        #print(o)
        offsets_relative(band, 4)
    else:
        log.setLevel('INFO')
        offsets_relative(band, 8)
