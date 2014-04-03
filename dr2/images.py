#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Prepares IPHAS images for public release.

This script will edit all images in the data release to ensure they have
the correct astrometric solution (WCS) and calibration information (PHOTZP)
in the header, which follows the conventions from the "ESO External Data
Products Standard" where possible
(cf. http://www.eso.org/sci/observing/phase3/p3edpstd.pdf )

Finally, the images are converted from Rice tile-compression into GZIP format,
which is slightly larger and slower but somewhat less obscure and in line with
many other archives.
"""
from __future__ import division, print_function, unicode_literals
from astropy.io import fits
from astropy.io import ascii
from astropy import log
from astropy import wcs
from astropy import table
import numpy as np
import itertools
import datetime
import os

import util
import constants

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Geert Barentsen']


####################
# CONSTANTS & CONFIG
####################

# Table containing slight updates to WCS astrometric parameters
WCSFIXES_PATH = os.path.join(constants.PACKAGEDIR, 'wcs-tuning', 'wcs-fixes.csv')
WCSFIXES = ascii.read(WCSFIXES_PATH)

MD = fits.getdata(os.path.join(constants.DESTINATION, 'metadata.fits'))
METADATA = dict(zip(MD['run'], MD))


###########
# CLASSES
###########

class CalibrationDatabase(object):
    """Class to hold the calibration shifts."""
    def __init__(self):
        r = ascii.read(os.path.join(constants.DESTINATION, 'calibration',
                                    'calibration-r.csv'))
        i = ascii.read(os.path.join(constants.DESTINATION, 'calibration',
                                    'calibration-i.csv'))
        ha = ascii.read(os.path.join(constants.DESTINATION, 'calibration',
                                     'calibration-ha.csv'))
        self.shifts = dict(zip(
                           np.concatenate((r['run'], i['run'], ha['run'])),
                           np.concatenate((r['shift'], i['shift'], ha['shift']))
                           ))


class SurveyImage(object):
    """Class used to write a single IPHAS CCD image with up-to-date keywords."""

    def __init__(self, run, ccd):
        self.run = run
        self.ccd = ccd
        # Open the image
        self.path_orig = constants.RAWDATADIR + METADATA[run]['image']
        self.fits_orig = fits.open(self.path_orig, do_not_scale_image_data=True)

        # Is the run a DR2-recalibrated run?
        mycaldb = get_caldb()
        if self.run in mycaldb.shifts:
            self.calibrated = True
        else:
            self.calibrated = False

        # Sort out the new FITS image and header
        self.hdu = fits.PrimaryHDU(self.fits_orig[self.ccd].data)
        self.set_header()
        self.fix_wcs()
        self.add_comments()

    def orig_header(self, keyword, extension=0):
        """Returns a keyword value from the original header."""
        try:
            return self.fits_orig[extension].header[keyword]
        except KeyError:
            return ""

    def orig_comments(self, keyword, extension=0):
        """Returns keyword comments from the original header."""
        try:
            return self.fits_orig[extension].header.comments[keyword]
        except KeyError:
            return ""

    @property
    def exptime(self):
        """Returns the exposure time as used in DR2.

        This can differ slightly from the original exposure time recorded in
        the raw data, because the DR2 pipeline accounts for known foibles
        in the in exposure time recording.
        """
        return METADATA[self.run]['exptime_precalib']

    @property
    def photzp(self):
        """Returns the zeropoint such that MAG=-2.5*log(pixel value)+PHOTZP

        Following the definition of PHOTZP defined in the
        "ESO External Data Products Standard"
        """
        # What is the calibration shift applied in DR2?
        mycaldb = get_caldb()
        try:
            shift = mycaldb.shifts[self.run]
        except KeyError:
            shift = 0.0
        # The zeropoint in the metadata file is corrected for extinction
        # but not re-calibrated and not corrected for PERCORR.
        # In accordance with the ESO standard, photzp absorbs the scaling
        # with exposure time.
        photzp = (METADATA[self.run]['zeropoint_precalib']
                  - self.percorr
                  + shift
                  + 2.5*np.log10(self.exptime))
        return np.round(photzp, 4)

    @property
    def percorr(self):
        return METADATA[self.run]['CCD{0}_PERCORR'.format(self.ccd)]

    @property
    def confmap(self):
        confmap = METADATA[self.run]['conf']
        if confmap == 'nan':
            return ''
        else:
            return confmap[1:]  # get rid of leading slash

    def set_header(self):
        # Copy keywords from the original HDU[0]
        for kw in ['RUN', 'OBSERVAT', 'OBSERVER', 'OBJECT',
                   'LATITUDE', 'LONGITUD', 'HEIGHT', 'SLATEL',
                   'TELESCOP',
                   'MJD-OBS', 'JD', 'PLATESCA', 'TELFOCUS',
                   'AIRMASS', 'DATE-OBS', 'UTSTART',
                   'TEMPTUBE', 'INSTRUME', 'WFFPOS', 'WFFBAND', 'WFFID',
                   'SECPPIX', 'DETECTOR', 'CCDSPEED',
                   'CCDXBIN', 'CCDYBIN', 'CCDSUM', 'CCDTEMP', 'NWINDOWS']:
            self.hdu.header[kw] = self.orig_header(kw)
            self.hdu.header.comments[kw] = self.orig_comments(kw)

        # Copy keywords from the original image extension
        for kw in ['BSCALE', 'BZERO', 'CCDNAME', 'CCDXPIXE', 'CCDYPIXE',
                   'AMPNAME', 'GAIN', 'READNOIS',
                   'NUMBRMS', 'STDCRMS',
                   'PERCORR', 'EXTINCT']:
            self.hdu.header[kw] = self.orig_header(kw, self.ccd)
            self.hdu.header.comments[kw] = self.orig_comments(kw, self.ccd)

        # Make it a proper ISO stamp
        self.hdu.header['DATE-OBS'] = self.fits_orig[0].header['DATE-OBS']+'T'+self.fits_orig[0].header['UTSTART']

        # Fix WCS - see fix_wcs() in detections.py!
        # Enforce the pipeline's defaults
        self.hdu.header['RADESYS'] = 'ICRS'
        self.hdu.header['EQUINOX'] = 2000.0
        self.hdu.header['CTYPE1'] = 'RA---ZPN'
        self.hdu.header['CTYPE2'] = 'DEC--ZPN'
        self.hdu.header['CRPIX1'] = self.fits_orig[self.ccd].header['CRPIX1']
        self.hdu.header['CRPIX2'] = self.fits_orig[self.ccd].header['CRPIX2']
        self.hdu.header['CRVAL1'] = self.fits_orig[self.ccd].header['CRVAL1']
        self.hdu.header['CRVAL2'] = self.fits_orig[self.ccd].header['CRVAL2']
        self.hdu.header['CRUNIT1'] = 'deg'
        self.hdu.header['CRUNIT2'] = 'deg'
        self.hdu.header['CD1_1'] = self.fits_orig[self.ccd].header['CD1_1']
        self.hdu.header['CD1_2'] = self.fits_orig[self.ccd].header['CD1_2']
        self.hdu.header['CD2_1'] = self.fits_orig[self.ccd].header['CD2_1']
        self.hdu.header['CD2_2'] = self.fits_orig[self.ccd].header['CD2_2']
        self.hdu.header['PV2_1'] = 1.0
        self.hdu.header['PV2_2'] = 0.0
        self.hdu.header['PV2_3'] = 220.0

        # Following the documentation at
        # http://apm49.ast.cam.ac.uk/surveys-projects/wfcam/technical/astrometry
        self.hdu.header.comments['RADESYS'] = 'WCS calibrated against 2MASS'
        self.hdu.header.comments['CTYPE1'] = 'Algorithm type for axis 1'
        self.hdu.header.comments['CTYPE2'] = 'Algorithm type for axis 2'
        self.hdu.header.comments['CRPIX1'] = '[pixel] Reference pixel along axis 1'
        self.hdu.header.comments['CRPIX2'] = '[pixel] Reference pixel along axis 2'
        self.hdu.header.comments['CRVAL1'] = '[deg] Right ascension at the reference pixel'
        self.hdu.header.comments['CRVAL2'] = '[deg] Declination at the reference pixel'
        self.hdu.header.comments['CRUNIT1'] = 'Unit of right ascension coordinates'
        self.hdu.header.comments['CRUNIT2'] = 'Unit of declination coordinates'
        self.hdu.header.comments['CD1_1'] = 'Transformation matrix element'
        self.hdu.header.comments['CD1_2'] = 'Transformation matrix element'
        self.hdu.header.comments['CD2_1'] = 'Transformation matrix element'
        self.hdu.header.comments['CD2_2'] = 'Transformation matrix element'
        self.hdu.header.comments['PV2_1'] = 'Coefficient for r term'
        self.hdu.header.comments['PV2_2'] = 'Coefficient for r**2 term'
        self.hdu.header.comments['PV2_3'] = 'Coefficient for r**3 term'

        # Fix zeropoint
        self.hdu.header['MAGZPT'] = self.orig_header('MAGZPT', self.ccd)
        self.hdu.header.comments['MAGZPT'] = 'Uncorrected nightly ZP (per second)'

        # Fix exposure time -- it might have changed in detections.py
        self.hdu.header['EXPTIME'] = self.exptime
        self.hdu.header.comments['EXPTIME'] = '[sec] Exposure time adopted in DR2'

        # True zeropoint with all corrections absorbed (including exposure time)
        self.hdu.header['PHOTZP'] = self.photzp
        self.hdu.header.comments['PHOTZP'] = 'mag(Vega) = -2.5*log(pixel value) + PHOTZP'

        # Add keywords according to the "ESO External Data Products standard"
        self.hdu.header['PHOTZPER'] = 0.03
        self.hdu.header.comments['PHOTZPER'] = 'Default 1-sigma PHOTZP uncertainty in IPHAS DR2'
        self.hdu.header['PHOTSYS'] = 'Vega'
        self.hdu.header.comments['PHOTSYS'] = 'Photometric system'

        # Was this image part of the DR2 re-calibration?
        if self.calibrated:
            self.hdu.header['FLUXCAL'] = 'ABSOLUTE'
        else:
            self.hdu.header['FLUXCAL'] = 'UNCALIBRATED'
        self.hdu.header.comments['FLUXCAL'] = 'Certifies the validity of PHOTZP'

        # Where is the confidence map?
        self.hdu.header['CONFMAP'] = self.confmap


    def fix_wcs(self):
        # Is an updated (fixed) WCS available?
        if self.run in WCSFIXES['RUN']:
            idx = ((WCSFIXES['RUN'] == self.run)
                   & (WCSFIXES['CCD'] == self.ccd))
            if idx.sum() > 0:
                log.info("WCS fixed: {0}[{1}].".format(self.run, self.ccd))
                idx_fix = idx.nonzero()[0][-1]
                for kw in ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                           'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
                    self.hdu.header[kw] = WCSFIXES[kw][idx_fix]


    def add_comments(self):
        """Populate the HISTORY and COMMENT keywords of the FITS file."""
        # Add history
        for line in str(self.fits_orig[self.ccd].header['HISTORY']).split('\n'):
            self.hdu.header['HISTORY'] = line
        self.hdu.header['HISTORY'] = datetime.datetime.now().strftime('%Y%m%d %H:%M:%S')
        self.hdu.header['HISTORY'] = 'Headers updated by Geert Barentsen as part of DR2.'
        self.hdu.header['HISTORY'] = 'This included changes to PHOTZP, EXPTIME and the WCS.'

        # Set calibration comments
        self.hdu.header['COMMENT'] = 'Photometric calibration info'
        self.hdu.header['COMMENT'] = '============================'

        if self.calibrated:
            self.hdu.header['COMMENT'] = 'The PHOTZP keyword in this header includes all the required'
            self.hdu.header['COMMENT'] = 'corrections for atmospheric extinction, gain variations,'
            self.hdu.header['COMMENT'] = 'exposure time, and the DR2 re-calibration shifts.'
            self.hdu.header['COMMENT'] = 'Hence to obtain calibrated magnitudes in the Vega system, use:'
            self.hdu.header['COMMENT'] = '    mag(Vega) = -2.5*log10(pixel value) + PHOTZP'
            self.hdu.header['COMMENT'] = 'If the brightness of a star is measured using a small aperture, then'
            self.hdu.header['COMMENT'] = 'you still need to add an aperture correction to this equation.'
        else:  # Uncalibrated image
            self.hdu.header['COMMENT'] = 'WARNING: this data is not part of DR2 and has not been re-calibrated.'
            self.hdu.header['COMMENT'] = 'It was likely excluded from DR2 due to a serious quality problem.'
            self.hdu.header['COMMENT'] = 'For example, PHOTZP might be inaccurate due to clouds.'
            self.hdu.header['COMMENT'] = '***USE AT YOUR OWN RISK***'

    @property
    def output_filename(self):
        """Filename of the output?"""
        return 'r{0}-{1}.fits.gz'.format(self.run, self.ccd).encode('ascii')

    def save(self):
        """Save the ccd image to a new file."""
        directory = os.path.join(constants.PATH_IMAGES,
                                 'r'+str(self.run)[0:3])
        util.setup_dir(directory)
        target = os.path.join(directory, self.output_filename)
        # checksum=True will add the CHECKSUM and DATASUM keywords
        self.hdu.writeto(target, clobber=True, checksum=True)

    def get_metadata(self):
        """Returns the CCD's metadata as a dictionary."""
        # Find center and corner coordinates (ra/dec in decimal degrees)
        mywcs = wcs.WCS(self.hdu.header)
        ra, dec = mywcs.all_pix2world([[1024, 2048]], 1)[0]
        corners = mywcs.all_pix2world([[1, 1],
                                       [1, 4096],
                                       [2048, 4096],
                                       [2048, 1]],
                                      1)
        ra1, ra2 = np.min(corners[:, 0]), np.max(corners[:, 0])
        dec1, dec2 = np.min(corners[:, 1]), np.max(corners[:, 1])

        if self.calibrated:
            in_dr2 = "true".encode('ascii')
        else:
            in_dr2 = "false".encode('ascii')

        band = str(self.hdu.header['WFFBAND']).lower()
        meta = {'filename': self.output_filename,
                'run': self.run,
                'ccd': self.ccd,
                'ra': ra,
                'dec': dec,
                'ra_min': ra1,
                'ra_max': ra2,
                'dec_min': dec1,
                'dec_max': dec2,
                'band': band,
                'utstart': str(self.hdu.header['DATE-OBS']).encode('ascii'),
                'exptime': self.hdu.header['EXPTIME'],
                'in_dr2': in_dr2,
                'seeing': METADATA[self.run]['CCD{0}_SEEING'.format(self.ccd)],
                'elliptic': METADATA[self.run]['CCD{0}_ELLIPTIC'.format(self.ccd)],
                'skylevel': METADATA[self.run]['CCD{0}_SKYLEVEL'.format(self.ccd)],
                'skynoise': METADATA[self.run]['CCD{0}_SKYNOISE'.format(self.ccd)],
                'airmass': self.hdu.header['AIRMASS'],
                'photzp': self.hdu.header['PHOTZP'],
                'confmap': self.confmap,
                }
        return meta


###########
# FUNCTIONS
###########

def get_caldb():
    """Returns the calibration information."""
    # Keep the CALDB stored as a global variable (= optimisation)
    global CALDB
    try:
        return CALDB
    except NameError:
        CALDB = CalibrationDatabase()
        return CALDB


def prepare_one(run):
    with log.log_to_file(os.path.join(constants.LOGDIR, 'images.log')):
        result = []
        for ccd in constants.EXTENSIONS:
            try:
                img = SurveyImage(run, ccd)
                img.save()
                result.append(img.get_metadata())
                img.fits_orig.close()  # avoid memory leak
            except Exception, e:
                log.error(str(run)+': '+util.get_pid()+': '+str(e))
        return result


def prepare_images(clusterview):
    # Make sure the output directory exists
    util.setup_dir(constants.PATH_IMAGES)
    metadata = []
    for band in ['halpha', 'r', 'i']:
        log.info('Starting with band {0}'.format(band))
        # Retrieve the list of runs
        if band == 'halpha':
            idx_band = 'ha'
        else:
            idx_band = band
        # [constants.IPHASQC_COND_RELEASE]
        runs = constants.IPHASQC['run_'+idx_band]
        # Prepare each run
        result = clusterview.map(prepare_one, runs[-600:-590], block=True)
        metadata.extend(result)

    # Write the metadata to a table
    mycolumns = (str('filename'), str('run'), str('ccd'),
                 str('ra'), str('dec'),
                 str('ra_min'), str('ra_max'),
                 str('dec_min'), str('dec_max'),
                 str('band'),
                 str('utstart'), str('exptime'),
                 str('in_dr2'),
                 str('seeing'), str('elliptic'),
                 str('skylevel'), str('skynoise'),
                 str('airmass'), str('photzp'),
                 str('confmap'))
    rows = list(itertools.chain.from_iterable(metadata))  # flatten list
    t = table.Table(rows, names=mycolumns)
    table_filename = os.path.join(constants.PATH_IMAGES, 'iphas-exposures.fits')
    t.write(table_filename, format='fits', overwrite=True)


################################
# MAIN EXECUTION (FOR DEBUGGING)
################################

if __name__ == '__main__':

    from IPython.parallel import client
    client = client.client.Client()
    with client[:].sync_imports():
        from dr2.images import SurveyImage
        from dr2 import constants
        from dr2 import util
        from astropy import log
        from astropy.io import fits
        import os
    prepare_images(client[:])

    #prepare_one(571408)
    #prepare_one(649055)
