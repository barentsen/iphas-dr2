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
from multiprocessing import Pool
import numpy as np
from astropy import log
from astropy.io import fits
import util
import constants
from constants import IPHASQC

np.seterr(invalid='ignore', divide='ignore')


class OffsetMachine(object):
    """Computes photometric offsets between exposures."""

    def __init__(self, run):
        self.run = run
        self.fits = fits.getdata(self.filename(run), 1)
        self.band = self.fits['band'][0]


    def filename(self, run):
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
        idx2 = (constants.COND_DR2
                & (dist < constants.FIELD_MAXDIST)
                & (IPHASQC['run_'+str(self.band)] != self.run))

        return IPHASQC['run_'+str(self.band)][idx2]

    def relative_offsets(self):
        offsets = []
        for run2 in self.overlap_runs():
            result = self._compute_relative_offsets(run2)
            offsets.append(result)
        return offsets

    def _compute_relative_offsets(self, run2):
        offsets = []
        fits2 = fits.getdata(self.filename(run2), 1)
        cond_reliable1 = ((self.fits['aperMag2'] > 15)
                          & (self.fits['aperMag2'] < 18)
                          & (self.fits['errBits'] == 0))
        cond_reliable2 = ((fits2['aperMag2'] > 14.8)
                          & (fits2['aperMag2'] < 18.2)
                          & (fits2['errBits'] == 0))

        for idx1 in np.argwhere(cond_reliable1):

            idx2 = crossmatch(self.fits['ra'][idx1],
                              self.fits['dec'][idx1],
                              fits2['ra'][cond_reliable2],
                              fits2['dec'][cond_reliable2])
            if idx2 is not None:
                offsets.append(fits2['aperMag2'][cond_reliable2][idx2]
                               - self.fits['aperMag2'][idx1])

        if len(offsets) < 5:
            return None
        else:
            offsets = np.array(offsets)
            return {'run1': self.run,
                    'run2': run2,
                    'offset': np.median(offsets),
                    'std': np.std(offsets),
                    'n': len(offsets)}


def crossmatch(ra, dec, ra_array, dec_array,
               matchdist=constants.MATCHING_DISTANCE):
    """Returns the index of the matched source."""
    dist = util.sphere_dist(ra, dec, ra_array, dec_array)
    idx_closest = dist.argmin()
    if dist[idx_closest] < (matchdist / 3600.):
        return idx_closest
    else:
        return None


def offsets_relative_one(run):
    try:
        om = OffsetMachine(run)
        offsets = om.relative_offsets()
        return offsets
    except Exception, e:
        log.error('UNEXPECTED EXCEPTION FOR RUN {0}: {1}'.format(run, e))
        return [None]

def offsets_relative(ncores=4):
    """Writes the file offsets-relative.csv with overlap offsets."""
    # Write the results
    filename = os.path.join(constants.DESTINATION, 'offsets-relative.csv')
    out = open(filename, 'w')
    out.write('run1,run2,offset,std,n\n')

    # Distribute the work over ncores
    p = Pool(processes=ncores)
    runs = np.concatenate((IPHASQC['run_r'][constants.COND_DR2],
                           IPHASQC['run_i'][constants.COND_DR2],
                           IPHASQC['run_ha'][constants.COND_DR2]))
    results = p.imap(offsets_relative_one, runs)
    for i, field in enumerate(results):
        log.info('Completed run {0}/{1}'.format(i, len(runs)))
        for row in field:
            if row is not None:
                out.write('{run1},{run2},{offset},{std},{n}\n'.format(**row))
        out.flush()

    out.close()


if __name__ == '__main__':

    log.setLevel('INFO')
    if constants.DEBUGMODE:
        #offsets_one(379750)
        offsets_relative(8)
    else:
        offsets_relative(8)

