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
from scipy import sparse
from scipy.sparse import linalg
from astropy.io import ascii
from astropy.io import fits
from astropy import log
import util
import constants
from constants import IPHASQC
from constants import IPHASQC_COND_RELEASE

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Geert Barentsen', 'Hywel Farnhill', 'Janet Drew']


PATH_UNCALIBRATED = os.path.join(constants.DESTINATION,
                                 'bandmerged')
PATH_CALIBRATED = os.path.join(constants.DESTINATION,
                               'bandmerged-calibrated')


# Extra anchors selected in the final phases of the data release,
# when a few areas with poor anchor coverage were spotted
EXTRA_ANCHORS = ['4510_jul2004a', '4510o_jul2004a',
                 '4583_jul2009', '4583o_jul2009',
                 '4525_jul2004a', '4525o_jul2004a',
                 '4598_jul2004a', '4598o_jul2004a',
                 '4087_jun2004',
                 '4099_jun2004', '4099o_jun2004',
                 '4110_jun2004', '4110o_jun2004',
                 '4096_jun2004', '4096o_jun2004']


class Calibration(object):
    """Holds the calibration shifts for all fields in the survey."""

    def __init__(self, band):
        """runs -- array of run identifiers"""
        #self.calib = np.array(zip(runs, np.zeros(len(runs))),
        #                      dtype=[('runs', 'i4'), ('shifts', 'f4')])
        assert(band in constants.BANDS)
        self.band = band

        self.runs = IPHASQC['run_'+band][IPHASQC_COND_RELEASE]
        self.shifts = np.zeros(len(self.runs))  # Shifts to be *ADDED*

        # APASS comparison data
        if band in ['r', 'i']:
            self.apass_shifts = IPHASQC[band+'shift_apassdr7'][IPHASQC_COND_RELEASE]
            self.apass_matches = IPHASQC[band+'match_apassdr7'][IPHASQC_COND_RELEASE]
        else:
            self.apass_shifts = np.zeros(len(self.runs))
            self.apass_matches = np.zeros(len(self.runs))

        assert(len(self.runs) == len(self.shifts))
        assert(len(self.runs) == len(self.apass_shifts))
        assert(len(self.apass_shifts) == len(self.apass_matches))

        self._load_offsetdata()

    def _load_offsetdata(self):
        filename_offsets = os.path.join(constants.DESTINATION,
                                        'offsets-{0}.csv'.format(self.band))
        log.info('Reading {0}'.format(filename_offsets))
        self.offsetdata = ascii.read(filename_offsets)

    def add_shifts(self, shifts):
        self.shifts += shifts

    def get_shift(self, run):
        """Shifts to be *ADDED* to the magnitudes of a run."""
        return self.shifts[self.runs == run][0]

    def evaluate(self):
        """Against APASS."""
        c_use = (self.apass_matches > 20)
        delta = self.apass_shifts[c_use] - self.shifts[c_use]
        stats =  "mean={0:.3f}+/-{1:.3f}, ".format(np.mean(delta),
                                                          np.std(delta))
        stats += "min/max={0:.3f}/{1:.3f}".format(np.min(delta),
                                                        np.max(delta))
        log.info(stats)

    def write(self, filename):
        """Write shifts to disk.
        """
        log.info('Writing results to {0}'.format(filename))
        f = open(filename, 'w')
        f.write('run,shift\n')
        for myrun, myshift in zip(self.runs, self.shifts):
            f.write('{0},{1}\n'.format(myrun, myshift))
        f.close()

    def get_overlaps(self):
        """Returns a dict with the magnitude offsets between run overlaps.

        Takes the current calibration into account.
        """
        log.info('Loading calibration-corrected magnitude offsets')
        # Performance optimisation
        current_shifts = dict(zip(self.runs, self.shifts))

        # Dictionary of field overlaps
        overlaps = {}
        for row in self.offsetdata:
            myrun1 = row['run1']
            myrun2 = row['run2']
            if myrun1 in self.runs and myrun2 in self.runs:
                # Offset is (run1 - run2), hence correcting for calibration
                # means adding (shift_run1 - shift_run2)
                myoffset = (row['offset'] 
                            + current_shifts[myrun1]
                            - current_shifts[myrun2])
                if myrun1 not in overlaps:
                    overlaps[myrun1] = {'runs': [], 'offsets': []}
                overlaps[myrun1]['runs'].append(myrun2)
                overlaps[myrun1]['offsets'].append(myoffset)

        return overlaps



def select_anchors():
    """Returns a boolean array indicating anchor status."""
    # 'anchors' is a boolean array indicating anchor status
    tolerance = 0.03  # Default = 0.03
    min_matches = 20  # Default = 20
    anchors = []
    APASS_OK = ( (IPHASQC.field('rmatch_apassdr7') >= min_matches)
                 & (IPHASQC.field('imatch_apassdr7') >= min_matches)
                 & (np.abs(IPHASQC.field('rshift_apassdr7')) <= tolerance)
                 & (np.abs(IPHASQC.field('ishift_apassdr7')) <= tolerance)
                 & (np.abs(IPHASQC.field('rshift_apassdr7') - IPHASQC.field('ishift_apassdr7')) <= tolerance) )

    APASS_ISNAN = ( np.isnan(IPHASQC.field('rshift_apassdr7'))
                    | np.isnan(IPHASQC.field('rshift_apassdr7'))
                    | (IPHASQC.field('rmatch_apassdr7') < min_matches)
                    | (IPHASQC.field('imatch_apassdr7') < min_matches) )

    SDSS_OK = ( (IPHASQC.field('rmatch_sdss') >= min_matches)
                 & (IPHASQC.field('imatch_sdss') >= min_matches)
                 & (np.abs(IPHASQC.field('rshift_sdss')) <= tolerance)
                 & (np.abs(IPHASQC.field('ishift_sdss')) <= tolerance)
                 & (np.abs(IPHASQC.field('rshift_sdss') - IPHASQC.field('ishift_sdss')) <= tolerance) )

    SDSS_ISNAN = ( np.isnan(IPHASQC.field('rshift_sdss'))
                    | np.isnan(IPHASQC.field('rshift_sdss'))
                    | (IPHASQC.field('rmatch_sdss') < min_matches)
                    | (IPHASQC.field('imatch_sdss') < min_matches) )

    # Allow extra anchors to be added
    EXTRA = np.array([myfield in EXTRA_ANCHORS for myfield in IPHASQC.field('id')])

    anchors = ( 
                ((IPHASQC.field('anchor') == 1) & APASS_ISNAN & (SDSS_OK | SDSS_ISNAN) )
                | APASS_OK
                | EXTRA
              )
    return anchors[IPHASQC_COND_RELEASE]

def select_anchors_halpha(runs):
    # Table containing information on photometric stability
    stabfile = os.path.join(constants.PACKAGEDIR,
                            'lib',
                            'photometric-stability.fits')
    stab = fits.getdata(stabfile, 1)

    anchors = []
    for i, run in enumerate(runs):
        # Is the run stable photometrically?
        cond_run = (stab['run_ha'] == run)
        if cond_run.sum() > 0:
            idx_run = np.argwhere(cond_run)[0][0]
            if ( (stab[idx_run]['hadiff'] > (-0.008-0.01))
                 and (stab[idx_run]['hadiff'] < (-0.008+0.01))
                 and (stab[idx_run]['rdiff'] > (-0.007-0.01))
                 and (stab[idx_run]['rdiff'] < (-0.007+0.01)) ):
                anchors.append(True)
            else:
                anchors.append(False)
    anchors = np.array(anchors)
    log.info('Found {0} H-alpha anchors'.format(anchors.sum()))
    return anchors


def calibrate_band(band='r'):
    """Calibrate a single band.

    band -- one of 'r', 'i', 'ha'
    """
    log.info('Starting to calibrate the {0} band'.format(band))

    """
    cal = Calibration(band)
    cal.evaluate()

    # Minimize overlap offsets
    anchors = select_anchors()
    overlaps = cal.get_overlaps()
    solver = Glazebrook(cal.runs, overlaps, anchors)    
    solver.solve()
    cal.add_shifts( solver.get_shifts() )
    cal.evaluate()
    plot_evaluation(cal, '/home/gb/tmp/cal-A.png', 'Strategy A')   
    """

    # H-alpha is a special case because the APASS-based selection of anchors
    # is not possible
    if band == 'ha':
        # We use the r-band calibration as the baseline for H-alpha
        rcalib_file = os.path.join(constants.DESTINATION,
                                   'calibration',
                                   'calibration-r.csv')
        rcalib = ascii.read(rcalib_file)
        cal = Calibration(band)
        cal.shifts = rcalib['shift']

        # We do run one iteration of Glazebrook using special H-alpha anchors
        overlaps = cal.get_overlaps()
        anchors = select_anchors_halpha(cal.runs)
        solver = Glazebrook(cal.runs, overlaps, anchors)    
        solver.solve()
        cal.add_shifts( solver.get_shifts() )

    else:
    
        cal = Calibration(band)
        plot_evaluation(cal, '{0}-step0.png'.format(band), 
                       '{0} - IPHAS-APASS uncalibrated'.format(band))

        # Correct outliers
        tolerance = 0.03  # Default = 0.03
        min_matches = 20  # Default = 20
        log.info('Correcting obvious outliers based on APASS')
        myshifts = np.zeros(len(cal.runs))
        c_outlier = ((cal.apass_matches > min_matches)
                     & (np.abs(cal.apass_shifts) > tolerance))
        myshifts[c_outlier] = cal.apass_shifts[c_outlier]
        cal.add_shifts(myshifts)
        cal.evaluate()
        plot_evaluation(cal, '{0}-step1.png'.format(band),
                        '{0} - Step 1: set initial conditions'.format(band))

        # Glazebrook: first pass (minimizes overlap offsets)
        anchors = select_anchors()
        overlaps = cal.get_overlaps()
        solver = Glazebrook(cal.runs, overlaps, anchors)    
        solver.solve()
        cal.add_shifts( solver.get_shifts() )
        cal.evaluate()
        plot_evaluation(cal, '{0}-step2.png'.format(band),
                       '{0} - Step 2: Glazebrook pass 1'.format(band))    
        
        # Add extra anchors
        delta = np.abs(cal.apass_shifts-cal.shifts)
        cond_extra_anchors = (cal.apass_matches > 20) & ~np.isnan(delta) & (delta > 0.05)
        log.info('Adding {0} extra anchors'.format(cond_extra_anchors.sum()))
        idx_extra_anchors = np.where(cond_extra_anchors)
        anchors = select_anchors()
        anchors[idx_extra_anchors] = True
        cal.shifts[idx_extra_anchors] = cal.apass_shifts[idx_extra_anchors]
        plot_evaluation(cal, '{0}-step3.png'.format(band),
                       '{0} - Step 3: added {1} extra anchors'.format(band, cond_extra_anchors.sum()))

        overlaps = cal.get_overlaps()
        solver = Glazebrook(cal.runs, overlaps, anchors)    
        solver.solve()
        cal.add_shifts( solver.get_shifts() )
        cal.evaluate()
        plot_evaluation(cal, '{0}-step4.png'.format(band),
                        '{0} - Step 4 - Glazebrook pass 2'.format(band))     
        

    filename = os.path.join(constants.DESTINATION, 
                            'calibration',
                            'calibration-{0}.csv'.format(band))
    cal.write(filename)

    return cal


def calibrate_test(band='r'):
    """Calibrate a single band.

    band -- one of 'r', 'i', 'ha'
    """
    log.info('Starting to calibrate the {0} band'.format(band))

    if band == 'ha':
        return Exception()

    else:
    
        cal = Calibration(band)

        # Minimize overlap offsets
        anchors = select_anchors()
        overlaps = cal.get_overlaps()
        solver = Glazebrook(cal.runs, overlaps, anchors)    
        solver.solve()
        cal.add_shifts( solver.get_shifts() )
        cal.evaluate()
        plot_evaluation(cal, '{0}-noinit.png'.format(band),
                       '{0} - Glazebrook without initial conditions changed'.format(band))    
    
    return cal


def plot_evaluation(cal,
                    filename, 
                    title='IPHAS-APASS after calibration'):
    c_use = (cal.apass_matches > 20)
    delta = cal.apass_shifts[c_use] - cal.shifts[c_use]
    l = IPHASQC['l'][IPHASQC_COND_RELEASE][c_use]
    b = IPHASQC['b'][IPHASQC_COND_RELEASE][c_use]

    import matplotlib
    matplotlib.use('Agg')  # Cluster does not have an X backend
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(12,6))
    fig.subplots_adjust(0.06, 0.15, 0.97, 0.9)
    p = fig.add_subplot(111)
    p.set_title(title)
    scat = p.scatter(l, b, c=delta, vmin=-0.13, vmax=+0.13,
                     edgecolors='none',
                     s=8, marker='h')
    plt.colorbar(scat)
    p.set_xlim([28, 217])
    p.set_ylim([-5.2, +5.2])
    p.set_xlabel('l')
    p.set_ylabel('b')

    path = os.path.join(constants.DESTINATION, 'calibration', filename)
    fig.savefig(path, dpi=200)
    log.info('Wrote {0}'.format(path))

    plt.close()
    return fig

def calibrate():
    # Make sure the output directory exists
    util.setup_dir(os.path.join(constants.DESTINATION, 'calibration'))
    
    for band in ['r', 'i', 'ha']:
        calibrate_band(band)



#############
# GLAZEBROOK
#############

class Glazebrook(object):
    """Finds zeropoints which minimise the offsets between overlapping fields.

    This class allows a set of catalogues with independently derived zeropoints
    to be brought to a global calibration with minimal magnitude offsets
    where fields overlap.

    This is achieved using the method detailed in the paper by 
    Glazebrook et al. 1994 (http://adsabs.harvard.edu/abs/1994MNRAS.266...65G).
    In brief, a set of equations are set up which allow the magnitude offsets
    between field overlaps to be minimised in a least squares sense.

    This class uses sparse matrix functions (scipy.sparse) to solve the large
    matrix equation in an efficient fashion.
    """

    def __init__(self, runs, overlaps, anchors):
        self.runs = runs
        self.overlaps = overlaps
        self.anchors = anchors
        self.nonanchors = ~anchors
        self.n_nonanchors = self.nonanchors.sum()
        log.info('There are {0} runs ({1} are anchors)'.format(len(runs),
                                                               anchors.sum()))

    def _A(self):
        """Returns the matrix called "A" in [Glazebrook 1994, Section 3.3]
        """
        log.info('Creating a sparse {0}x{0} matrix (might take a while)'.format(self.n_nonanchors))
        A = sparse.lil_matrix((self.n_nonanchors,
                               self.n_nonanchors))

        nonanchorruns = self.runs[self.nonanchors]
        # Loop over all non-anchors that make up the matrix
        for i, run in enumerate(nonanchorruns):
            overlapping_runs = self.overlaps[run]['runs']
            
            # On the diagonal, the matrix holds the negative number of overlaps
            A[i, i] = -float(len(overlapping_runs))

            # Off the diagonal, the matrix holds a 1 where two runs overlap
            for run2 in overlapping_runs:
                idx_run2 = np.argwhere(run2 == nonanchorruns)
                if len(idx_run2) > 0:
                    j = idx_run2[0]  # Index of the overlapping run
                    A[i, j] = 1.
                    A[j, i] = 1.     # Symmetric matrix
        return A

    def _b(self):
        """Returns the vector called "b" in [Glazebrook 1994, Section 3.3]
        """
        b = np.zeros(self.n_nonanchors)
        for i, run in enumerate(self.runs[self.nonanchors]):
            b[i] = np.sum(self.overlaps[run]['offsets'])
        return b

    def solve(self):
        """Returns the solution of the matrix equation.
        """
        self.A = self._A()
        self.b = self._b()
        log.info('Now solving the matrix equation')
        # Note: there may be alternative algorithms
        # which are faster for symmetric matrices.
        self.solution = linalg.lsqr(self.A, self.b,
                                    atol=1e-8, iter_lim=2e5, show=False)
        log.info('Solution found')
        log.info('mean shift = {0} +/- {1}'.format(
                                            np.mean(self.solution[0]),
                                            np.std(self.solution[0])))
        return self.solution

    def get_shifts(self):
        shifts = np.zeros(len(self.runs))
        shifts[self.nonanchors] = self.solution[0]
        return shifts

    def write(self, filename):
        """Write shifts to disk.
        """
        log.info('Writing results to {0}'.format(filename))
        f = open(filename, 'w')
        f.write('run,shift\n')
        for i, myrun in enumerate(self.runs[self.nonanchors]):
            f.write('{0},{1}\n'.format(myrun, self.solution[0][i]))
        for myrun in self.runs[self.anchors]:
            f.write('{0},{1}\n'.format(myrun, 0.0))
        f.close()


class CalibrationApplicator(object):
    """Applies the calibration to the catalogues.

    This class will read a bandmerged catalogue from the 'bandmerged' directory,
    apply the appropriate calibration shifts as listed in 
    'calibration/calibration-{r,i,ha}.csv', and then write the updated 
    catalogue to a new directory 'bandmerged-calibrated'."""

    def __init__(self):
        self.datadir = PATH_UNCALIBRATED
        self.outdir = PATH_CALIBRATED
        util.setup_dir(self.outdir)

        # Read in the calibration
        self.calib = {}
        for band in constants.BANDS:
            calib_file = os.path.join(constants.DESTINATION,
                                      'calibration',
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
                        == IPHASQC.field('run_' + band)[idx_field])
            if cond_run.sum() > 0:
                shifts[band] = self.calib[band]['shift'][cond_run][0]
            else:
                log.warning('No shift for %s' % filename)
                shifts[band] = 0.0

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


def calibrate_one(filename):
    """Applies the photometric re-calibration to a single bandmerged field catalogue."""
    with log.log_to_file(os.path.join(constants.LOGDIR, 'apply_calibration.log')):
        try:
            ca = CalibrationApplicator()
            ca.run(filename)
        except Exception, e:
            log.error('%s: *UNEXPECTED EXCEPTION*: calibrate_one: %s' % (filename, e))
        return filename


def apply_calibration(clusterview):
    """Applies the photometric re-calibration to all bandmerged field catalogues."""
    filenames = os.listdir(PATH_UNCALIBRATED)
    results = clusterview.imap(calibrate_one, filenames)

    # Print a friendly message once in a while
    i = 0
    for filename in results:
        i += 1
        if (i % 1000) == 0:
            log.info('Completed file {0}/{1}'.format(i, len(filenames)))
    log.info('Application of calibration finished')


def median_rmha_one(path):
    """Returns the median of a single band-merged catalogue."""
    mydata = fits.getdata(path, 1)
    mask_reliable = mydata['reliable']
    fieldid = path.split('.')[-2]
    median = np.median(mydata['rmha'][mask_reliable])
    return (fieldid, median)

def median_rmha(clusterview,
                directory=PATH_UNCALIBRATED,
                output_filename = os.path.join(constants.DESTINATION,
                                              'calibration',
                                              'median-rmha.csv')):
    """Computes the median r-Ha colour in all fields.

    This will be used as an input to evaluate the H-alpha calibration.
    """
    log.info('Starting to compute median(r-ha) values.')
    util.setup_dir(os.path.join(constants.DESTINATION, 'calibration')
    paths = [os.path.join(directory, filename) 
             for filename in os.listdir(directory)]
    results = clusterview.imap(median_rmha_one, paths[0:10])
    # Write the results to a csv file
    with open(output_filename, 'w') as out:
        out.write('field,median_rmha\n')
        for (myfield, mymedian) in results:
            out.write('{0},{1}\n'.format(myfield, mymedian))
    log.info('Computing median(r-ha) values finished.')    


###################
# MAIN EXECUTION
###################

if __name__ == '__main__':
    log.setLevel('DEBUG')
    calibrate()
