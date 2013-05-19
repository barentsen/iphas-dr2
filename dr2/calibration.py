#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fits a global photometric calibration using the Glazebrook algorithm.

The algorithm finds a set of zeropoint shifts which minimizes the magnitude
offsets between overlapping exposures (computed using the dr2.offsets module.)

This file also contains a class to apply the calibration to the catalogues.

TODO
* Calibrate on a CCD-by-CCD basis?
"""
import numpy as np
import os
from multiprocessing import Pool
from scipy import sparse
from scipy.sparse import linalg
from astropy.io import ascii
from astropy.io import fits
from astropy import log
import constants
from constants import IPHASQC

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Hywel Farnhill', 'Janet Drew']


#############
# GLAZEBROOK
#############

class Glazebrook(object):
    """Fit a global calibration using the method by Glazebrook et al. (1994)

    Glazebrook et al. (1994) is hereafter referred to as [Glazebrook]
    Uses sparse matrix functions (scipy.sparse) for efficiency.
    """

    def __init__(self, runs, overlaps, anchors):
        self.runs = runs
        self.overlaps = overlaps
        self.anchors = anchors
        log.info('There are {0} runs ({1} are anchors)'.format(len(runs),
                                                               anchors.sum()))
        self.n_nonanchors = (~anchors).sum()

    def _A(self):
        """Returns the matrix called "A" in [Glazebrook]
        """
        log.info('Creating a sparse {0}x{0} matrix'.format(self.n_nonanchors))
        A = sparse.lil_matrix((self.n_nonanchors,
                               self.n_nonanchors))
        for i, run in enumerate(self.runs[~self.anchors]):
            for j, run2 in enumerate(self.runs[~self.anchors]):
                if j < i:  # Symmetric matrix
                    continue
                elif i == j:
                    A[i, j] = -len(self.overlaps[run2]['runs'])
                elif run2 in self.overlaps[run]['runs']:
                    A[i, j] = 1
                    A[j, i] = 1  # Symmetric matrix
        return A

    def _b(self):
        """Returns the vector called "b" in [Glazebrook]
        """
        b = np.zeros(self.n_nonanchors)
        for i, run in enumerate(self.runs[~self.anchors]):
            b[i] = np.sum(self.overlaps[run]['offsets'])
        return b

    def solve(self):
        """Returns the solution of the matrix equation.
        """
        self.A = self._A()
        self.b = self._b()
        log.info('Now solving the matrix equation')
        # Note: there should be alternative algorithms for symmetric
        # matrices which are faster.
        self.solution = linalg.lsqr(self.A, self.b,
                                    atol=1e-10, iter_lim=2e5, show=False)
        #log.info('Solution found: {0}'.format(self.solution))
        log.info('Solution found')
        return self.solution

    def write(self, filename):
        """Write shifts to disk.
        """
        log.info('Writing results to {0}'.format(filename))
        f = open(filename, 'w')
        f.write('run,shift\n')
        for i, myrun in enumerate(self.runs[~self.anchors]):
            f.write('{0},{1}\n'.format(myrun, self.solution[0][i]))
        for myrun in self.runs[self.anchors]:
            f.write('{0},{1}\n'.format(myrun, 0.0))
        f.close()


def glazebrook_data(band='r'):
    """Retrieve the offset information needed by the Glazebrook algorithm.
    """
    assert(band in ['r', 'i', 'ha'])

    # Parse data
    filename_offsets = os.path.join(constants.DESTINATION,
                                    'offsets-{0}.csv'.format(band))
    log.info('Reading {0}'.format(filename_offsets))
    offsetdata = ascii.read(filename_offsets)

    # All the runs, sorted numerically
    runs = np.sort(np.unique(offsetdata['run1']))

    # 'anchors' is a boolean array indicating anchor status
    anchors = []

    cond_anchors = ( IPHASQC.field('is_best')
                     & ( (IPHASQC.field('anchor') == 1)
                         | ( (np.abs(IPHASQC.field('rshift_apassdr7')) <= 0.03)
                            & (np.abs(IPHASQC.field('ishift_apassdr7')) <= 0.03) )
                        )
                     )

    log.info('Found {0} anchors'.format(cond_anchors.sum()))


    """
    cond_anchors = (IPHASQC.field('anchor') == 1)
    """

    QC_RUNS = IPHASQC.field('run_{0}'.format(band))
    for myrun in runs:
        if cond_anchors[QC_RUNS == myrun][0]:
            anchors.append(True)
        else:
            anchors.append(False)
    anchors = np.array(anchors)

    # Dictionary of field overlaps
    overlaps = {}
    for row in offsetdata:
        myrun = row['run1']
        if myrun not in overlaps:
            overlaps[myrun] = {'runs': [], 'offsets': []}
        overlaps[myrun]['runs'].append( row['run2'] )
        overlaps[myrun]['offsets'].append( row['offset'] )

    return (runs, overlaps, anchors)


def run_glazebrook_band(band='r'):
    """Produce the calibration for a given band.
    """
    assert(band in ['r', 'i', 'ha'])
    log.info('Running Glazebrook for the {0} band'.format(band))
    filename_output = os.path.join(constants.DESTINATION,
                                   'calibration-{0}.csv'.format(band))
    runs, overlaps, anchors = glazebrook_data(band)
    g = Glazebrook(runs, overlaps, anchors)
    g.solve()
    g.write(filename_output)


def run_glazebrook(ncores=3):
    p = Pool(processes=ncores)
    p.map(run_glazebrook_band, constants.BANDS)


###################################################
# CLASS USED TO APPLY THE CALIBRATION TO CATALOGUES
###################################################

PATH_UNCALIBRATED = os.path.join(constants.DESTINATION,
                                 'bandmerged')
PATH_CALIBRATED = os.path.join(constants.DESTINATION,
                               'bandmerged-calibrated')


class CalibrationApplicator(object):
    """Updates the bandmerged catalogues by applying the calibration shifts."""

    def __init__(self):
        self.datadir = PATH_UNCALIBRATED
        self.outdir = PATH_CALIBRATED
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Read in the calibration
        self.calib = {}
        for band in constants.BANDS:
            calib_file = os.path.join(constants.DESTINATION,
                                      'calibration-{0}.csv'.format(band))
            self.calib[band] = ascii.read(calib_file)

    def run(self, filename):
        #for filename in os.listdir(self.datadir):
        log.info('Correcting {0}'.format(filename))
        self.calibrate(filename)

    def get_shifts(self, filename):
        fieldid = filename.split('.fits')[0]
        idx_field = np.argwhere(IPHASQC.field('id') == fieldid)[0]

        shifts = {}
        for band in constants.BANDS:
            cond_run = (self.calib[band]['run']
                        == IPHASQC.field('run_'+band)[idx_field])
            shifts[band] = self.calib[band]['shift'][cond_run][0]

        log.info("Shifts for {0}: {1}".format(fieldid, shifts))
        return shifts

    def calibrate(self, filename):
        path_in = os.path.join(self.datadir, filename)
        path_out = os.path.join(self.outdir, filename)
        shifts = self.get_shifts(filename)

        param = {'stilts': constants.STILTS,
                 'filename_in': path_in,
                 'filename_out': path_out,
                 'cmd': """'replacecol r "toFloat(r  + {r})"; \
                            replacecol rPeakMag "toFloat(rPeakMag  + {r})"; \
                            replacecol rAperMag1 "toFloat(rAperMag1  + {r})"; \
                            replacecol rAperMag3 "toFloat(rAperMag3  + {r})"; \
                            replacecol i "toFloat(i  + {i})"; \
                            replacecol iPeakMag "toFloat(iPeakMag  + {i})"; \
                            replacecol iAperMag1 "toFloat(iAperMag1  + {i})"; \
                            replacecol iAperMag3 "toFloat(iAperMag3  + {i})"; \
                            replacecol ha "toFloat(ha  + {ha})"; \
                            replacecol haPeakMag "toFloat(haPeakMag  + {ha})"; \
                            replacecol haAperMag1 "toFloat(haAperMag1  + {ha})"; \
                            replacecol haAperMag3 "toFloat(haAperMag3  + {ha})"; \
                            replacecol rmi "toFloat(r-i)"; \
                            replacecol rmha "toFloat(r-ha)"; \
                           '""".format(**shifts)}

        cmd = '{stilts} tpipe cmd={cmd} in={filename_in} out={filename_out}'.format(**param)
        log.debug(cmd)
        status = os.system(cmd)
        log.info('stilts status: '+str(status))
        return status


"""
def apply_calibration_strip(strip):
    log.info('Applying calibration to strip {0}'.format(strip))
    ca = CalibrationApplicator(strip)
    ca.run()

def apply_calibration(ncores=2):
    strips = np.arange(25, 215+1, constants.STRIPWIDTH)[::-1]
    p = Pool(processes=ncores)
    p.map(apply_calibration_strip, strips)
"""

def calibration_worker(filename):
    ca = CalibrationApplicator()
    ca.run(filename)


def apply_calibration(ncores=2):
    filenames = os.listdir(PATH_UNCALIBRATED)
    p = Pool(processes=ncores)
    p.map(calibration_worker, filenames)


def evaluate_calibration(band='r'):
    """Returns the array of residuals of the calibration against APASS"""
    # Validation data
    filename_apass = os.path.join(constants.PACKAGEDIR,
                                  'tests',
                                  'apass_validation.fits')
    apass = fits.getdata(filename_apass, 1)

    filename_calib = os.path.join(constants.DESTINATION,
                                  'calibration-{0}.csv'.format(band))
    calib = ascii.read(filename_calib)

    filename_eval = os.path.join(constants.DESTINATION,
                                 'eval-{0}.csv'.format(band))
    out = open(filename_eval, 'w')
    out.write('id,night,l,b,n_apass,apass,calib\n')
    residuals = []
    for i in range(len(calib)):
        idx = np.argwhere(apass['run_{0}'.format(band)] == calib['run'][i])[0]
        apass_shift = apass['{0}shift_apassdr7'.format(band)][idx][0]
        calib_shift = calib['shift'][i]

        apass_stars = apass['{0}match_apassdr7'.format(band)][idx][0]

        out.write('{0},{1},{2},{3},{4},{5},{6}\n'.format(apass['id'][idx][0],
                                                    apass['night'][idx][0],
                                                    apass['l'][idx][0],
                                                    apass['b'][idx][0],
                                                    apass_stars,
                                                    apass_shift,
                                                    calib_shift))
        if not np.isnan(apass_shift):
            residuals.append( apass_shift - calib_shift )

    out.close()

    print('===== {0} ====='.format(band))
    print('mean(residuals {0}): {1:.03f}'.format(band, np.mean(residuals)))
    print('min(residuals {0}): {1:.03f}'.format(band, np.min(residuals)))
    print('max(residuals {0}): {1:.03f}'.format(band, np.max(residuals)))
    print('std(residuals {0}): {1:.03f}'.format(band, np.std(residuals)))
    return residuals
    #crossmatch using stilts

###################
# MAIN EXECUTION
###################

if __name__ == '__main__':

    if constants.DEBUGMODE:
        log.setLevel('DEBUG')
        #run_glazebrook(ncores=3)
        #apply_calibration_strip(215)
        bands = ['r', 'i']
        p = Pool(processes=2)
        p.map(run_glazebrook_band, bands)
        for band in bands:
            residuals = evaluate_calibration(band)

    else:
        log.setLevel('INFO')
        run_glazebrook(ncores=3)
        #apply_calibration(ncores=8)
        apply_calibration(8)
