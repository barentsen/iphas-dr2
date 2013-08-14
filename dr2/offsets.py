#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Computes magnitude offsets between overlapping exposures.

This script carries out pair-wise crossmatching of reliable objects in
overlapping exposures. It then computes the median offset between those
magnitudes, which provides the input required to carry out a global
re-calibration.

Output
------
The output is a CSV file written to 
"{constants.DESTINATION}/calibration/offsets-{band}.csv"

The columns in the CSV file are:
run1   -- reference run number
run2   -- comparison run number
offset -- median(run1_magnitudes - run2_magnitudes)
std    -- std(run1_magnitudes - run2_magnitudes)
n      -- len(run1_magnitudes - run2_magnitudes)

Dependencies
------------
* IPHASQC table containing all metadata.
* single-band detection catalogues.
"""
from __future__ import division, print_function, unicode_literals
import os
import sys
from multiprocessing import Pool
import numpy as np
from astropy import log
from astropy.io import fits

import constants
from constants import IPHASQC
import util

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Geert Barentsen', 'Hywel Farnhill', 'Janet Drew']


#############################
# CONSTANTS & CONFIGURATION
#############################

# Magnitude ranges to consider reliable
MAGLIMITS = {'r': [15, 17.5], 'i': [14, 16.5], 'ha': [15, 17.5]}

# Maximum matching distance
MATCHING_DISTANCE = 0.5  # arcsec


###########
# CLASSES
###########

class OffsetMachine(object):
    """Computes photometric offsets between exposures."""

    def __init__(self, run):
        """
        Parameters
        ----------
        run : string or int
            Telescope exposure identifier.
        """
        self.run = run
        self.data = self.get_data(run)
        self.band = self.data['band']

    def get_data(self, myrun):
        """Returns a dictionary with the data for a given run.

        Parameters
        ----------
        myrun : string or int
            Telescope exposure identifier.

        Returns
        -------
        data : dictionary of arrays
        """
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
        """Returns the list of overlapping exposures in the same band.

        Returns
        -------
        runs : list
            List of overlapping exposure identifiers.
        """
        cond_run = (IPHASQC['run_'+str(self.band)] == self.run)
        idx = np.argwhere(cond_run)[0]
        myra = IPHASQC['ra'][idx]
        mydec = IPHASQC['dec'][idx]

        dist = util.sphere_dist(myra, mydec, IPHASQC['ra'], IPHASQC['dec'])
        idx2 = ( constants.IPHASQC_COND_RELEASE
                    & (dist < constants.FIELD_MAXDIST)
                    & (IPHASQC['run_'+str(self.band)] != self.run)
                )

        return IPHASQC['run_'+str(self.band)][idx2]

    def relative_offsets(self):
        """Returns the offsets for all runs that overlap with self.run.

        Returns
        -------
        offsets : list of dictionaries.
            Each dictionary contains the offset information for an overlapping 
            exposure.
        """
        log.debug('Computing offsets for run {0}'.format(self.run))
        offsets = []
        for run2 in self.overlap_runs():
            log.debug(str(run2))
            result = self._compute_relative_offsets(run2)
            if result is not None:
                offsets.append(result)
        return offsets

    def _compute_relative_offsets(self, run2):
        """Computes the offset against a specified run.

        Parameters
        ----------
        run2 : string or int
            The run to compare against.

        Returns
        -------
        offsets : dictionary
            Containing the fields run1, run2, offset, std, n.
        """
        limit_bright = MAGLIMITS[self.band][0]
        limit_faint = MAGLIMITS[self.band][1]

        offsets = []

        offset_data = self.get_data(run2)

        cond_reliable1 = ((self.data['aperMag2'] > limit_bright)
                          & (self.data['aperMag2'] < limit_faint)
                          & (self.data['errBits'] == 0))
        cond_reliable2 = ((offset_data['aperMag2'] > (limit_bright-0.2))
                          & (offset_data['aperMag2'] < (limit_faint+0.2))
                          & (offset_data['errBits'] == 0))

        for idx1 in np.argwhere(cond_reliable1):
            idx2 = util.crossmatch(self.data['ra'][idx1],
                                   self.data['dec'][idx1],
                                   offset_data['ra'][cond_reliable2],
                                   offset_data['dec'][cond_reliable2],
                                   matchdist=MATCHING_DISTANCE)
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

def offsets_one(run):
    """Returns all offsets for a given reference exposure.

    Parameters
    ----------
    run : integer or string
        Telescope run number.

    Returns
    -------
    offsets : list of dictionaries
        A sequence of dictionaries for each overlapping run. 
        Each dictionary contains the fields run1/run2/offset/std/n.
    """
    with log.log_to_file(os.path.join(constants.LOGDIR, 'offsets.log')):
        try:
            log.info('{0}: Computing offsets for {1}'.format(util.get_pid(), run)
            om = OffsetMachine(run)
            return om.relative_offsets()
        except Exception, e:
            log.error('{0}: UNEXPECTED EXCEPTION FOR RUN {1}: {2}'.format(util.get_pid(),
                                                                          run,
                                                                          e))
            return [None]


def compute_offsets_band(clusterview, band, 
                         destination=os.path.join(constants.DESTINATION,
                                                  'calibration')):
    """Computes magnitude offsets between all overlapping runs in a given band.

    The output is a file called offsets-{band}.csv which contains the columns
        run1   -- reference exposure (telescope run number)
        run2   -- comparison exposure (telescope run number)
        offset -- median(run1_magnitudes - run2_magnitudes)
        std    -- stdev(run1_magnitudes - run2_magnitudes)
        n      -- number of crossmatched stars used in computing offset/std.

    Parameters
    ----------
    clusterview : cluster view derived used e.g. IPython.parallel.Client()[:]
        Work will be spread across the nodes in this cluster view.

    band : string
        One of 'r', 'i', 'ha'.

    destination : string
        Directory where the output csv file will be written.
    """
    assert(band in constants.BANDS)
    log.info('Starting to compute offsets for band {0}'.format(band))

    # Write the results
    util.setup_dir(destination)
    filename = os.path.join(destination, 'offsets-{0}.csv'.format(band))
    out = open(filename, 'w')
    out.write('run1,run2,offset,std,n\n')

    # Distribute the work across the cluster
    runs = IPHASQC['run_'+str(band)][constants.IPHASQC_COND_RELEASE]
    #runs = IPHASQC['run_'+str(band)]
    results = clusterview.imap(offsets_one, runs)

    # Write offsets to the CSV file as the results are returned
    i = 0
    for offsets in results:
        i += 1
        for row in offsets:
            if row is not None:
                out.write('{run1},{run2},{offset},{std},{n}\n'.format(**row))
        # Print a friendly status message once in a while
        if (i % 100) == 0:
            log.info('Completed run {0}/{1}'.format(i, len(runs)))
            out.flush()

    out.close()


def compute_offsets(clusterview):
    """Computes magnitude offsets for all overlapping runs and all bands.

    Parameters
    ----------
    clusterview : cluster view derived used e.g. IPython.parallel.Client()[:]
        Work will be spread across the nodes in this cluster view.
    """   
    for band in constants.BANDS:
        compute_offsets_band(clusterview, band)

