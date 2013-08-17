#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fits a global photometric calibration using the Glazebrook algorithm.

The algorithm finds a set of zeropoint shifts which minimizes the magnitude
offsets between overlapping exposures (computed using the dr2.offsets module.)
In addition, the APASS survey is used to ensure that we do not deviate from
the 'absolute truth'.

This is the most complicated and tricky part of the data release procedure,
because IPHAS suffered from a fair number of nights with variable transparency
and hence zeropoint uncertainties.

This file also contains a class to apply the calibration to the catalogues,
and functions to plot colour-colour diagrams for quality control.
"""
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')  # Cluster does not have an X backend
from matplotlib import pyplot as plt
from scipy import sparse
from scipy.sparse import linalg
from astropy.io import ascii
from astropy.io import fits
from astropy import log

import constants
from constants import IPHASQC
from constants import IPHASQC_COND_RELEASE
from constants import CALIBDIR
import util

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Geert Barentsen', 'Hywel Farnhill', 'Janet Drew']


############
# CONSTANTS
############

# When to trust other surveys?
TOLERANCE = 0.03 # abs(iphas-apass) tolerated
MIN_MATCHES = 30 # minimum number of matches in a field against reference survey


##########
# CLASSES
##########

class Calibration(object):
    """Container for calibration information in a single band.

    This class holds information about the offsets between overlaps, 
    the offsets against other surveys, the choice of anchor fields, 
    and the calibration shifts required.

    This class is effectively a container to hold everything we know
    about our survey zeropoints, and contains several functions to 
    interact with this information (e.g. create spatial plots of zeropoint 
    offsets.)

    Attributes
    ----------
    band : {'r', 'i', 'ha'}
    runs : array of int
        List of exposure numbers for `band` which are part of the data release.
    shifts : array of float
        The calibration shifts to be *added* to the magnitudes of `runs`.
    anchors : array of bool
        Which exposures can be trusted?
    """

    def __init__(self, band):
        """Loads the necessary information about the survey zeropoints.

        Parameters
        ----------
        band : str {'r', 'i', 'ha'}
            Name of the photometric filter being calibrated.
        """
        #self.calib = np.array(zip(runs, np.zeros(len(runs))),
        #                      dtype=[('runs', 'i4'), ('shifts', 'f4')])
        assert(band in constants.BANDS)
        self.band = band

        self.runs = IPHASQC['run_'+band][IPHASQC_COND_RELEASE]
        self.shifts = np.zeros(len(self.runs))  # Shifts to be *ADDED* - init to 0

        # Load broad-band comparison data
        if band in ['r', 'i']:
            self.apass_shifts = IPHASQC[band+'shift_apassdr7'][IPHASQC_COND_RELEASE]
            self.apass_matches = IPHASQC[band+'match_apassdr7'][IPHASQC_COND_RELEASE]
            self.sdss_shifts = IPHASQC[band+'shift_sdss'][IPHASQC_COND_RELEASE]
            self.sdss_matches = IPHASQC[band+'match_sdss'][IPHASQC_COND_RELEASE]
        else:
            self.apass_shifts = np.zeros(len(self.runs))
            self.apass_matches = np.zeros(len(self.runs))

        # Sanity check: we should have the same number of runs and shifts 
        assert(len(self.runs) == len(self.shifts))
        assert(len(self.runs) == len(self.apass_shifts))
        assert(len(self.apass_shifts) == len(self.apass_matches))

        self._load_offsetdata()
        self.anchors = self.select_anchors()

    def get_runs(self):
        return self.runs

    def get_anchors(self):
        return self.anchors

    def get_overlaps(self, weights=True):
        """Returns a dict with the magnitude offsets between run overlaps.

        Takes the current calibration into account.
        """
        log.info('Loading calibration-corrected magnitude offsets between overlaps')
        # Prepare look-up dictionary which links runs to calibration shifts
        current_shifts = dict(zip(self.runs, self.shifts))

        # Dictionary of field overlaps
        overlaps = {}
        for row in self.offsetdata:
            myrun1 = row['run1']
            myrun2 = row['run2']
            if myrun1 in self.runs and myrun2 in self.runs:
                # Offset is computed as (run1 - run2), hence correcting 
                # for calibration means adding (shift_run1 - shift_run2)
                myoffset = (row['offset'] 
                            + current_shifts[myrun1]
                            - current_shifts[myrun2])
                if myrun1 not in overlaps:
                    overlaps[myrun1] = {'runs': [], 'offsets': [], 'weights': []}
                overlaps[myrun1]['runs'].append(myrun2)
                overlaps[myrun1]['offsets'].append(myoffset)

                if weights:
                    overlaps[myrun1]['weights'].append(np.sqrt(row['n']))
                else:
                    overlaps[myrun1]['weights'].append(1.0)

        return overlaps

    def write_anchor_list(self, filename):
        """Writes the list of anchors to a csv files.

        Parameters
        ----------
        filename : str
        """
        with open(filename, 'w') as out:
            out.write('run,is_anchor\n')
            for i in range(len(self.runs)):
                out.write('{0},{1}\n'.format(self.runs[i], self.anchors[i]))

    def add_shifts(self, shifts):
        self.shifts += shifts

    def get_shift(self, run):
        """Returns the calibrations shift for a given run.

        Parameters
        ----------
        run : integer
            Exosure identifier for which you want to know the calibration shift.

        Returns
        -------
        shift : float
            Shift to be *added* to the magnitudes of the specified run.
        """
        return self.shifts[self.runs == run][0]

    def evaluate(self, name, title):
        # Plot the absolute calibration shifts
        l = IPHASQC['l'][IPHASQC_COND_RELEASE]
        b = IPHASQC['b'][IPHASQC_COND_RELEASE]
        self._spatial_plot(l, b, self.shifts, 'calib-'+name, 'Calibration '+title)

        if self.band in ['r', 'i']:
            statsfile = os.path.join(CALIBDIR, 'stats-{0}.txt'.format(self.band))
            with open(statsfile, 'w') as out:
                # Against APASS
                mask_use = (self.apass_matches >= MIN_MATCHES)
                l = IPHASQC['l'][IPHASQC_COND_RELEASE][mask_use]
                b = IPHASQC['b'][IPHASQC_COND_RELEASE][mask_use]
                delta = self.apass_shifts[mask_use] - self.shifts[mask_use]
                self._spatial_plot(l, b, delta, 'apass-'+name, 'APASS: '+title)

                stats =  "mean={0:.3f}+/-{1:.3f}, ".format(np.mean(delta),
                                                           np.std(delta))
                stats += "min/max={0:.3f}/{1:.3f}".format(np.min(delta),
                                                          np.max(delta))

                out.write(stats)
                log.info(stats)

                # Against SDSS
                mask_use = (self.sdss_matches >= MIN_MATCHES)
                l = IPHASQC['l'][IPHASQC_COND_RELEASE][mask_use]
                b = IPHASQC['b'][IPHASQC_COND_RELEASE][mask_use]
                delta = self.sdss_shifts[mask_use] - self.shifts[mask_use]
                self._spatial_plot(l, b, delta, 'sdss-'+name, 'SDSS '+title)


    def _spatial_plot(self, l, b, shifts, name, title=''):
        """Creates a spatial plot of l/b against shifts."""
        plotdir = os.path.join(CALIBDIR, 'plots')
        util.setup_dir(plotdir)

        fig = plt.figure(figsize=(12,6))
        fig.subplots_adjust(0.06, 0.15, 0.97, 0.9)
        p = fig.add_subplot(111)
        p.set_title(title)
        scat = p.scatter(l, b, c=shifts, vmin=-0.13, vmax=+0.13,
                         edgecolors='none',
                         s=7, marker='h')
        plt.colorbar(scat)
        p.set_xlim([28, 217])
        p.set_ylim([-5.2, +5.2])
        p.set_xlabel('l')
        p.set_ylabel('b')

        path = os.path.join(plotdir, self.band+'-'+name+'-without-anchors.png')
        fig.savefig(path, dpi=200)
        log.info('Wrote {0}'.format(path))

        # Indicate anchors
        p.scatter(IPHASQC['l'][IPHASQC_COND_RELEASE][self.anchors],
                  IPHASQC['b'][IPHASQC_COND_RELEASE][self.anchors],
                  edgecolors='black', facecolor='none',
                  s=15, marker='x', alpha=0.9, lw=0.3)
        path = os.path.join(plotdir, self.band+'-'+name+'-with-anchors.png')
        fig.savefig(path, dpi=200)
        log.info('Wrote {0}'.format(path))

        plt.close()
        return fig

    def write(self, filename):
        """Writes calibration shifts to a CSV file on disk.

        Parameters
        ----------
        filename : string
            Filename of the CSV file to write the calibration shifts.
        """
        log.info('Writing results to {0}'.format(filename))
        f = open(filename, 'w')
        f.write('run,shift\n')
        for myrun, myshift in zip(self.runs, self.shifts):
            f.write('{0},{1}\n'.format(myrun, myshift))
        f.close()

    def _load_offsetdata(self):
        filename_offsets = os.path.join(CALIBDIR,
                                        'offsets-{0}.csv'.format(self.band))
        log.info('Reading {0}'.format(filename_offsets))
        mydata = ascii.read(filename_offsets)
        # Do not use the offsets unless enough stars were used
        #mask_use = (mydata['n'] >= 5) & (mydata['std'] < 0.1)
        self.offsetdata = mydata

    def select_anchors(self):
        """Returns a boolean array indicating which runs are suitable anchors."""
        median_pair_offset = -0.008 
        IS_STABLE = ( (IPHASQC.field('med_dr') < (median_pair_offset+0.03)) &
                      (IPHASQC.field('med_dr') > (median_pair_offset-0.03)) &
                      (IPHASQC.field('med_di') < (median_pair_offset+0.03)) &
                      (IPHASQC.field('med_di') > (median_pair_offset-0.03)) &
                      (IPHASQC.field('med_dh') < (median_pair_offset+0.03)) &
                      (IPHASQC.field('med_dh') > (median_pair_offset-0.03))
                    )

        if self.band == 'ha':
            # Because the H-alpha calibration is tied to the r-band,
            # we require fields to be "stable" to be an anchor in H-alpha.
            # "Stable" is defined as the fieldpair not showing great shifts.
            anchors = IS_STABLE[IPHASQC_COND_RELEASE]
            log.info('IS_STABLE: {0} fields are H-alpha anchors'.format(anchors.sum()))
            return anchors

        else:
            tolerance = 0.03  # Default = 0.03
            min_matches = 30  # Default = 20
            anchors = []
            IS_APASS_ANCHOR = ( (IPHASQC.field('rmatch_apassdr7') >= min_matches)
                              & (IPHASQC.field('imatch_apassdr7') >= min_matches)
                              & (np.abs(IPHASQC.field('rshift_apassdr7')) <= tolerance)
                              & (np.abs(IPHASQC.field('ishift_apassdr7')) <= tolerance)
                              & (np.abs(IPHASQC.field('rshift_apassdr7') - IPHASQC.field('ishift_apassdr7')) <= tolerance) )

            IS_OLD_ANCHOR = (IPHASQC.field('anchor') == 1)

            # Extra anchors selected in the final phases of the data release,
            # when a few areas with poor anchor coverage were spotted
            EXTRA_ANCHORS = ascii.read('lib/anchor-extra.txt')['field']
            IS_EXTRA_ANCHOR = np.array([myfield in EXTRA_ANCHORS 
                                        for myfield in IPHASQC.field('id')])
            # Make sure the following runs are no anchors
            # cf. e-mail Janet to Geert, 13 Aug 2013
            ANCHOR_BLACKLIST = ascii.read('lib/anchor-blacklist.txt')['field']
            IS_BLACKLIST = np.array([myfield in ANCHOR_BLACKLIST 
                                     for myfield in IPHASQC.field('id')])

            # Extra night with good conditions
            ANCHOR_NIGHTS = [20030915, 20031018, 20031101, 20031104, 20031108,
                             20031117, 20040707, 20040805, 20040822, 20041022,
                             20050629, 20050709, 20050710, 20050711, 20050916,
                             20050917, 20050918, 20051023, 20051101, 20051102,
                             20061129, 20061130, 20061214, 20070627, 20070630,
                             20080722, 20080723, 20090808, 20090810, 20091029,
                             20091031]
            IS_IN_EXTRA_NIGHT = np.array([mynight in ANCHOR_NIGHTS 
                                       for mynight in IPHASQC.field('night')])

            # Nights which should NOT provide anchors
            NIGHT_BLACKLIST = [20031117, 20051109, 20061128, 20091029, 20101029,]
            IS_IN_NIGHT_BLACKLIST = np.array([night in NIGHT_BLACKLIST 
                                              for night in IPHASQC.field('night')])

            # Anchors must not have known quality issues
            IS_QUALITY_OK = ((IPHASQC.field('seeing_max') < 2.0) &
                               (IPHASQC.field('airmass_max') < 1.4) &
                               (IPHASQC.field('qflag') != 'C') &
                               (IPHASQC.field('qflag') != 'D'))

            anchors = (-IS_BLACKLIST &
                       -IS_IN_NIGHT_BLACKLIST &
                       IS_STABLE & 
                       IS_QUALITY_OK &
                       (IS_OLD_ANCHOR | IS_EXTRA_ANCHOR | IS_IN_EXTRA_NIGHT | IS_APASS_ANCHOR)
                      )
            result = anchors[IPHASQC_COND_RELEASE]

            log.info('IS_APASS_ANCHOR: {0} fields'.format(IS_APASS_ANCHOR.sum()))
            log.info('IS_OLD_ANCHOR: {0} fields'.format(IS_OLD_ANCHOR.sum()))
            log.info('IS_EXTRA_ANCHOR: {0} fields'.format(IS_EXTRA_ANCHOR.sum()))
            log.info('IS_IN_EXTRA_NIGHT: {0} fields'.format(IS_IN_EXTRA_NIGHT.sum()))
            log.info('IS_BLACKLIST: {0} fields'.format(IS_BLACKLIST.sum()))            
            log.info('Anchors in data release: {0} fields'.format(result.sum()))

            return result


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

    def __init__(self, cal):
        self.cal = cal  # Calibration object
        self.runs = cal.get_runs()
        self.overlaps = cal.get_overlaps()
        self.anchors = cal.get_anchors()

        self.nonanchors = ~self.anchors
        self.n_nonanchors = self.nonanchors.sum()
        log.info('Glazebrook: there are {0} runs ({1} are anchors)'.format(
                                                            len(self.runs),
                                                            self.anchors.sum()))

    def _A(self):
        """Returns the matrix called "A" in [Glazebrook 1994, Section 3.3]
        """
        log.info('Glazebrook: creating a sparse {0}x{0} matrix (might take a while)'.format(self.n_nonanchors))
        A = sparse.lil_matrix((self.n_nonanchors,
                               self.n_nonanchors))

        nonanchorruns = self.runs[self.nonanchors]
        # Loop over all non-anchors that make up the matrix
        for i, run in enumerate(nonanchorruns):
            
            try:
                # On the diagonal, the matrix holds the negative sum of weights
                A[i, i] = -float(np.sum(self.overlaps[run]['weights']))

                # Off the diagonal, the matrix holds the weight where two runs overlap
                for run2, weight in zip(self.overlaps[run]['runs'],
                                        self.overlaps[run]['weights']):
                    idx_run2 = np.argwhere(run2 == nonanchorruns)
                    if len(idx_run2) > 0:
                        j = idx_run2[0]  # Index of the overlapping run
                        A[i, j] = weight
                        A[j, i] = weight  # Symmetric matrix
            except KeyError:
                log.warning('Glazebrook: no overlap data for run {0}'.format(run))
                A[i, i] = -1.0
        return A

    def _b(self):
        """Returns the vector called "b" in [Glazebrook 1994, Section 3.3]
        """
        b = np.zeros(self.n_nonanchors)
        for i, run in enumerate(self.runs[self.nonanchors]):
            try:
                b[i] = np.sum(
                              np.array(self.overlaps[run]['offsets']) *
                              np.array(self.overlaps[run]['weights'])
                              )
            except KeyError:
                log.warning('Glazebrook: no overlap data for run {0}'.format(run))
                b[1] = 0.0
        return b

    def solve(self):
        """Returns the solution of the matrix equation.
        """
        self.A = self._A()
        self.b = self._b()
        log.info('Glazebrook: now solving the matrix equation')
        # Note: there may be alternative algorithms
        # which are faster for symmetric matrices.
        self.solution = linalg.lsqr(self.A, self.b,
                                    atol=1e-8, iter_lim=1e6, show=False)
        log.info('Glazebrook: solution found')
        log.info('Glazebrook: mean shift = {0} +/- {1}'.format(
                                            np.mean(self.solution[0]),
                                            np.std(self.solution[0])))

        shifts = np.zeros(len(self.runs))
        shifts[self.nonanchors] = self.solution[0]
        return shifts


class CalibrationApplicator(object):
    """Applies the calibration to a bandmerged catalogue.

    This class will read a bandmerged catalogue from the 'bandmerged' 
    directory, apply the appropriate calibration shifts as listed in 
    'calibration/calibration-{r,i,ha}.csv', and then write the updated 
    catalogue to a new directory called 'bandmerged-calibrated'.
    """

    def __init__(self):
        self.datadir = constants.PATH_BANDMERGED
        self.outdir = constants.PATH_BANDMERGED_CALIBRATED
        util.setup_dir(self.outdir)

        # Read the calibration information into a dictionary
        self.calib = {}
        for band in constants.BANDS:
            calib_file = os.path.join(CALIBDIR, 'calibration-{0}.csv'.format(band))
            self.calib[band] = ascii.read(calib_file)

    def run(self, filename):
        #for filename in os.listdir(self.datadir):
        log.info('Correcting {0}'.format(filename))
        self.calibrate(filename)

    def get_shifts(self, filename):
        """Retuns the calibration shifts for a given bandmerged catalogue.
        """
        fieldid = filename.split('.fits')[0]
        return get_field_shifts(fieldid)


    def get_field_shifts(self, fieldid):
        """Returns the calibration shifts for a given field.

        Parameters
        ----------
        fieldid : str
            Field identifier, e.g. "0001_aug2003"

        Returns
        -------
        shifts : dictionary {'r':shift_r, 'i': shift_i, 'ha': shift_ha}
            Shifts to add to the magnitudes to calibrate a field.
        """
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

        log.debug("Shifts for {0}: {1}".format(fieldid, shifts))
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


##
## Apply the calibration to the bandmerged catalogues
##

def apply_calibration(clusterview):
    """Applies the photometric re-calibration to all bandmerged field catalogues."""
    filenames = os.listdir(constants.PATH_BANDMERGED)
    results = clusterview.imap(calibrate_one, filenames)

    # Print a friendly message once in a while
    i = 0
    for filename in results:
        i += 1
        if (i % 1000) == 0:
            log.info('Completed file {0}/{1}'.format(i, len(filenames)))
    log.info('Application of calibration finished')

def calibrate_one(filename):
    """Applies the photometric re-calibration to a single bandmerged field catalogue."""
    with log.log_to_file(os.path.join(constants.LOGDIR, 'apply_calibration.log')):
        try:
            ca = CalibrationApplicator()
            ca.run(filename)
        except Exception, e:
            log.error('%s: *UNEXPECTED EXCEPTION*: calibrate_one: %s' % (filename, e))
        return filename


##
## Plotting colour-colour diagrams of anchors and fields for quality control
##

def plot_anchors():
    """Plots diagrams of the anchors."""
    # Setup output directory
    inputdir = constants.PATH_BANDMERGED
    outputdir = os.path.join(CALIBDIR, 'anchors')
    util.setup_dir(outputdir)
    # Which fields to plot?
    anchorlist = ascii.read(os.path.join(CALIBDIR, 'anchors-r-initial.csv'))
    anchor_runs = anchorlist['run'][anchorlist['is_anchor'] == 'True']
    fields = [util.run2field(myrun, 'r') for myrun in anchor_runs]
    # Distribute work
    log.info('Starting to plot {0} anchors'.format(len(fields)))
    from multiprocessing import Pool
    mypool = Pool(4)
    mypool.map(plot_field, zip(fields,
                               [inputdir]*len(fields),
                               [outputdir]*len(fields),
                               [{'r':0.0, 'i':0.0, 'ha':0.0}]*len(fields)))


"""
def plot_calibrated_fields():
    # Plots diagrams of all fields after calibration.
    # Setup output directory
    inputdir = constants.PATH_BANDMERGED_CALIBRATED
    outputdir = os.path.join(CALIBDIR, 'diagrams')
    util.setup_dir(outputdir)
    # Which fields to plot?
    fields = IPHASQC['id'][IPHASQC_COND_RELEASE]
    # Distribute work
    from multiprocessing import Pool
    mypool = Pool(4)
    mypool.map(plot_field, zip(fields,
                               [inputdir]*len(fields),
                               [outputdir]*len(fields)))
"""

def plot_calibrated_fields():
    inputdir = constants.PATH_BANDMERGED
    outputdir = os.path.join(CALIBDIR, 'diagrams')
    util.setup_dir(outputdir)
    ca = CalibrationApplicator()
    args = []
    for i, field in enumerate(IPHASQC['id'][IPHASQC_COND_RELEASE]):
        args.append((field, inputdir, outputdir, ca.get_field_shifts(field)))
    
    log.info('Starting to plot {0} anchors'.format(len(args)))
    from multiprocessing import Pool
    mypool = Pool(4)
    mypool.map(plot_field, args)


def plot_field(arguments):
    field, inputdir, outputdir, shifts = arguments
    shift_rmi = shifts['r'] - shifts['i']
    shift_rmha = shifts['r'] - shifts['ha']
    # Initiate figure
    fig = plt.figure(figsize=(6,4))
    fig.subplots_adjust(0.15, 0.15, 0.95, 0.9)
    p = fig.add_subplot(111)
    p.set_title('{0} (r {1:+.2f}, i {2:+.2f}, ha {3:+.2f})'.format(
                field, shifts['r'], shifts['i'], shifts['ha']),
                fontsize=14)
    # Load colours from the bandmerged catalogue
    d = fits.getdata(os.path.join(inputdir, field+'.fits'), 1)
    mask_use = (d['r'] < 19.0) & (d['errBits'] == 0) & (d['pStar'] > 0.2)
    p.scatter(d['rmi'][mask_use] + shift_rmi,
              d['rmha'][mask_use] + shift_rmha,
              alpha=0.4, edgecolor="red", facecolor="red",
              lw=0, s=1, marker='o')
    # Main sequence
    p.plot([0.029, 0.212, 0.368, 0.445, 0.903, 1.829],
           [0.001, 0.114, 0.204, 0.278, 0.499, 0.889],
           c='black', lw=0.5)
    # A-type reddening line
    p.plot([0.029, 0.699, 1.352, 1.991, 2.616],
           [0.001, 0.199, 0.355, 0.468, 0.544],
           c='black', lw=0.5)

    p.set_xlim([-0.2, +2.0])
    p.set_ylim([-0.1, +1.3])
    p.set_xlabel('r-i')
    p.set_ylabel('r-Ha')
    # Write to disk
    path = os.path.join(outputdir, field+'.jpg')
    fig.savefig(path, dpi=200)
    plt.close()
    log.info('Wrote {0}'.format(path))




##
## The functions which drive the zeropoint calibration
##

def calibrate_band(band='r'):
    """Calibrate a single band.

    Parameters
    ----------
    band : one of 'r', 'i', 'ha'

    Returns
    -------
    cal : Calibration class
        object entailing the shifts to be added (cal.shifts)
    """
    log.info('Starting to calibrate the {0} band'.format(band))

    # H-alpha is a special case because the APASS-based selection of anchors
    # is not possible
    if band == 'ha':
        # We use the r-band calibration as the baseline for H-alpha
        rcalib = ascii.read(os.path.join(CALIBDIR, 'calibration-r.csv'))
        cal = Calibration(band)
        cal.shifts = rcalib['shift']
        cal.evaluate('step1', 'H-alpha with r-band shifts')

        # We do run one iteration of Glazebrook using special H-alpha anchors
        overlaps = cal.get_overlaps()
        solver = Glazebrook(cal)
        shifts = solver.solve()
        cal.add_shifts(shifts)
        cal.evaluate('step2', 'H-alpha after Glazebrook')

        cal.write_anchor_list(os.path.join(CALIBDIR, 'anchors-{0}-initial.csv'.format(band)))

    else:
    
        cal = Calibration(band)


        # Hack: take account of exposure time changes
        MANUALLY_SHIFTED = [364687,
            368903, 368904, 368923, 368925,
            369998, 370073, 370076, 370084,
            370095, 371652, 371695, 372557,
            372684, 372707, 372751, 372771,
            372880, 373106, 373111, 373698,
            374904, 376449, 376461, 376463,
            376481, 376493, 376530, 401548,
            401566, 402270, 407505, 407580,
            407586, 407598, 408287, 408296,
            413548, 413566, 413596, 413783,
            413804, 414671, 418169, 418190,
            418196, 418310, 427588, 427820,
            457662, 460468, 470277, 470592,
            470822, 470852, 474652, 476050,
            476131, 478320, 478434, 478609,
            478645, 478720, 478795, 537478,
            537544, 537550, 537565, 537623,
            538318, 538354, 538366, 538406,
            538595, 538601, 538759, 540932,
            541185, 541717, 541948, 568871,
            568892, 568937, 568970, 568982,
            569666, 569768, 569816, 
            570005, 570559, 570601, 570754,
            571311, 571362, 571377, 571704,
            597412, 597469, 597778, 598536,
            598710, 598865, 598880, 647562,
            649761, 686153, 686264, 687199,
            687757, 702703, 702724, 702769,
            703360, 703408, 703741]
        MANUALLY_SHIFTED = np.array([myrun in SHIFTED_RUNS for myrun in cal.runs])

        cal.evaluate('step1', '{0} - uncalibrated'.format(band))
        cal.write_anchor_list(os.path.join(CALIBDIR, 'anchors-{0}-initial.csv'.format(band)))

        # Glazebrook: first pass (minimizes overlap offsets)
        solver = Glazebrook(cal)
        shifts = solver.solve()
        cal.add_shifts(shifts)
        cal.evaluate('step2', '{0} - step 2: Glazebrook pass 1'.format(band))    
        # Write the used anchors to a csv file


        # Correct outliers against APASS and fix them as anchors
        delta = np.abs(cal.apass_shifts - cal.shifts)
        cond_extra_anchors = ((cal.apass_matches >= MIN_MATCHES) &
                              -np.isnan(delta) &
                              (delta >= TOLERANCE) &
                              -MANUALLY_SHIFTED)
        idx_extra_anchors = np.where(cond_extra_anchors)
        cal.anchors[idx_extra_anchors] = True
        cal.shifts[idx_extra_anchors] = cal.apass_shifts[idx_extra_anchors]
        log.info('Adding {0} extra anchors'.format(cond_extra_anchors.sum()))
        cal.evaluate('step3',
                     '{0} - step 3: added {1} extra anchors'.format(
                                       band, cond_extra_anchors.sum()))

        # Run Glazebrook again with the newly added anchors
        solver = Glazebrook(cal)    
        shifts = solver.solve()
        cal.add_shifts(shifts)
        cal.evaluate('step4', '{0} - step 4 - Glazebrook pass 2'.format(band))

        # Remove bad APASS shifts one more time
        mask_has_apass = ((cal.apass_matches >= MIN_MATCHES) &
                          -np.isnan(delta))
        mask_is_outlier = (mask_has_apass &
                           (np.abs(cal.apass_shifts - cal.shifts) >= TOLERANCE) &
                           -MANUALLY_SHIFTED
                           )
        idx_outlier = np.where(mask_is_outlier)
        cal.shifts[idx_outlier] = cal.apass_shifts[idx_outlier]
        # And make all APASS overlaps anchors
        cal.anchors[mask_has_apass] = True 
        cal.evaluate('step5', '{0} - step 5: anchor all APASS fields'.format(band))

        # Run Glazebrook again with the newly added anchors
        solver = Glazebrook(cal)    
        shifts = solver.solve()
        cal.add_shifts( shifts )
        cal.evaluate('step6', '{0} - step 6 - Glazebrook pass 3'.format(band))
        

    # This is the important file: the calibration shifts!
    filename = os.path.join(CALIBDIR, 'calibration-{0}.csv'.format(band))
    cal.write(filename)

    return cal


def calibrate():
    """Calibrates all bands in the survey.

    Produces files called "calibration{r,i,ha}.csv" which tabulate
    the zeropoint shifts to be *added* to each exposure.
    """
    # Make sure the output directory exists
    util.setup_dir(CALIBDIR)
    # Calibrate each band in the survey
    for band in constants.BANDS:
        calibrate_band(band)


def calibrate_multiprocessing():
    """Calibrates all bands in the survey.

    Produces files called "calibration{r,i,ha}.csv" which tabulate
    the zeropoint shifts to be *added* to each exposure.
    """
    # Make sure the output directory exists
    util.setup_dir(CALIBDIR)
    # Calibrate each band in the survey
    from multiprocessing import Pool
    pool = Pool(2)
    pool.map(calibrate_band, ['r', 'i'])
    calibrate_band('ha') # H-alpha depends on output of r


################################
# MAIN EXECUTION (FOR DEBUGGING)
#################################

if __name__ == '__main__':
    log.setLevel('DEBUG')
    calibrate()
    #calibrate_multiprocessing()
    #calibrate_band('ha')
    plot_anchors()
    plot_calibrated_fields()
