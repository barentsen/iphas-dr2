#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calibrates IPHAS detection catalogues against APASS/SDSS.

Tasks:
- Find shift against APASS & SDSS (average).
- Find shifts amongst fields.
- Solve Glazebrook algorithm in 5x5 strips.
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
__credits__ = ['Hywel Farnhill', 'Robert Greimel', 'Janet Drew']


#############################
# CONSTANTS & CONFIGURATION
#############################

# Magnitude limits to compute offsets from
MAGLIMITS = {'r': [15, 17.5], 'i': [14, 16.5], 'ha': [15, 17.5]}


###########
# CLASSES
###########

class OffsetMachine(object):
    """Computes photometric offsets between exposures."""

    def __init__(self, run):
        self.run = run
        self.fits = fits.getdata(self.filename(run), 1)
        self.band = self.fits['band'][0]

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
        idx2 = (constants.COND_RELEASE
                & (dist < constants.FIELD_MAXDIST)
                & (IPHASQC['run_'+str(self.band)] != self.run))

        return IPHASQC['run_'+str(self.band)][idx2]

    def relative_offsets(self):
        log.debug('Computing offsets for run {0}'.format(self.run))
        offsets = []
        for run2 in self.overlap_runs():
            result = self._compute_relative_offsets(run2)
            if result is not None:
                offsets.append(result)
        return offsets

    def _compute_relative_offsets(self, run2):
        limit_bright = MAGLIMITS[self.band][0]
        limit_faint = MAGLIMITS[self.band][1]

        offsets = []
        fits2 = fits.getdata(self.filename(run2), 1)
        cond_reliable1 = ((self.fits['aperMag2'] > limit_bright)
                          & (self.fits['aperMag2'] < limit_faint)
                          & (self.fits['errBits'] == 0))
        cond_reliable2 = ((fits2['aperMag2'] > limit_bright-0.2)
                          & (fits2['aperMag2'] < limit_faint+0.2)
                          & (fits2['errBits'] == 0))

        for idx1 in np.argwhere(cond_reliable1):

            idx2 = util.crossmatch(self.fits['ra'][idx1],
                                   self.fits['dec'][idx1],
                                   fits2['ra'][cond_reliable2],
                                   fits2['dec'][cond_reliable2])
            if idx2 is not None:
                offsets.append(self.fits['aperMag2'][idx1]
                               - fits2['aperMag2'][cond_reliable2][idx2])

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
        offsets = om.relative_offsets()
        return offsets
    except Exception, e:
        log.error('UNEXPECTED EXCEPTION FOR RUN {0}: {1}'.format(run, e))
        return [None]

def offsets_relative(band, ncores=4):
    """Writes the file offsets-relative.csv with overlap offsets."""
    assert(band in constants.BANDS)
    log.info('Starting to compute offsets for band {0}'.format(band))

    # Write the results
    filename = os.path.join(constants.DESTINATION, 'offsets-r.csv')
    out = open(filename, 'w')
    out.write('run1,run2,offset,std,n\n')

    # Distribute the work over ncores
    p = Pool(processes=ncores)
    runs = IPHASQC['run_'+str(band)][constants.COND_RELEASE]
    results = p.map(offsets_relative_one, runs[0:4])
    for i, field in enumerate(results):
        log.info('Completed run {0}/{1}'.format(i, len(runs)))
        for row in field:
            if row is not None:
                out.write('{run1},{run2},{offset},{std},{n}\n'.format(**row))

        if (i % 50) == 0:  # Make sure we write at least every 50 runs
            out.flush()

    out.close()


if __name__ == '__main__':

    # Which band to process?
    if len(sys.argv) > 1:
        band = sys.argv[1]
    else:
        raise Exception('Please give the band as the first argument')

    if constants.DEBUGMODE:
        log.setLevel('DEBUG')
        o = offsets_relative_one(571029)
        print(o)
        #offsets_relative(band, 4)
    else:
        offsets_relative(band, 8)
