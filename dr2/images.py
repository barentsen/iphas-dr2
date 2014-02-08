#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Prepare IPHAS images for data release.

This script will edit all images in the data release to ensure they have
the correct astrometric solution (WCS) and calibration information (zeropoint)
in the header.

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

# Where the write output images?
IMAGE_DESTINATION = os.path.join(constants.DESTINATION, 'images')

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

CALDB = CalibrationDatabase()


class SurveyImage(object):
    """Class used to write a single IPHAS CCD image with up-to-date keywords."""

    def __init__(self, run, ccd):
        self.run = run
        self.ccd = ccd
        # Open the image
        self.path_orig = constants.RAWDATADIR + METADATA[run]['image']
        self.fits_orig = fits.open(self.path_orig, do_not_scale_image_data=True)

        # Is the run a DR2-recalibrated run?
        if self.run in CALDB.shifts:
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
    def output_filename(self):
        """Filename of the output?"""
        return '{0}_{1}.fits.gz'.format(self.run, self.ccd).encode('ascii')

    @property
    def exptime(self):
        """Returns the exposure time as used in DR2.

        This can differ slightly from the original exposure time recorded in
        the raw data, because the DR2 pipeline accounts for known foibles
        in the in exposure time recording.
        """
        return METADATA[self.run]['exptime_precalib']

    @property
    def zeropoint(self):
        """Returns the image's zeropoint ZP such that `mag = ZP - 2.5*log(flux)`
        """
        # What is the calibration shift applied in DR2?
        try:
            shift = CALDB.shifts[self.run]
        except KeyError:
            shift = 0.0
        # The zeropoint in the metadata file is corrected for extinction
        # but not re-calibrated and not corrected for PERCORR
        return METADATA[self.run]['zeropoint_precalib'] - self.percorr + shift

    @property
    def percorr(self):
        return METADATA[self.run]['CCD{0}_PERCORR'.format(self.ccd)]

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
        self.hdu.header['ORIGZPT'] = self.orig_header('MAGZPT', self.ccd)
        self.hdu.header.comments['ORIGZPT'] = 'Original nightly ZP; uncorrected for extinction/clouds'
        
        self.hdu.header['MAGZPT'] = self.zeropoint
        if self.calibrated:
            self.hdu.header.comments['MAGZPT'] = 'Re-calibrated DR2 zeropoint'
        else:
            self.hdu.header.comments['MAGZPT'] = 'ORIGZPT corrected for exinction and PERCORR'

        # Fix exposure time -- it might have changed in detections.py
        self.hdu.header['EXPTIME'] = self.exptime
        self.hdu.header.comments['EXPTIME'] = '[sec] Exposure time assumed by the pipeline'

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
        self.hdu.header['HISTORY'] = '    Headers updated by Geert Barentsen as part of DR2.'
        self.hdu.header['HISTORY'] = '    This included changes to MAGZPT, EXPTIME and the WCS.'

        # Set calibration comments
        self.hdu.header['COMMENT'] = 'Calibration info'
        self.hdu.header['COMMENT'] = '================'

        if self.calibrated:
            self.hdu.header['COMMENT'] = 'The MAGZPT keyword in this header has been corrected for atmospheric'
            self.hdu.header['COMMENT'] = 'extinction and gain (PERCORR) and has been re-calibrated as part of DR2.'
            self.hdu.header['COMMENT'] = ''
            self.hdu.header['COMMENT'] = 'Hence to obtain calibrated magnitudes relative to Vega, use:'
            self.hdu.header['COMMENT'] = '    mag(Vega) = MAGZPT - 2.5*log(pixel value / EXPTIME)'
        else:
            self.hdu.header['COMMENT'] = 'Warning: this data is not part of DR2 and has not been re-calibrated.'
            self.hdu.header['COMMENT'] = 'It was likely excluded from DR2 for a serious quality problem.'
            self.hdu.header['COMMENT'] = 'Use at your own risk.'

    def save(self):
        """Save the ccd image to a new file."""
        directory = os.path.join(IMAGE_DESTINATION,
                                 str(self.hdu.header['WFFBAND']).lower())
        target = os.path.join(directory, self.output_filename)
        self.hdu.writeto(target, clobber=True)

    def get_metadata(self):
        """Returns the CCD's metadata as a dictionary."""
        # Find center and corner coordinates (ra/dec in decimal degrees)
        mywcs = wcs.WCS(self.hdu.header)
        ra, dec = mywcs.all_pix2world([[1024, 2048]], 1)[0]
        corners = mywcs.all_pix2world([[1,1], [1,4096], [2048,4096], [2048,1]], 1)
        ra1, ra2 = np.min(corners[:,0]), np.max(corners[:,0])
        dec1, dec2 = np.min(corners[:,1]), np.max(corners[:,1])
        
        if self.calibrated:
            in_dr2 = "true".encode('ascii')
        else:
            in_dr2 = "false".encode('ascii')

        meta = {'filename': self.output_filename,
                'band': str(self.hdu.header['WFFBAND']).lower(),
                'dr2': in_dr2,
                'run': self.run,
                'ccd': self.ccd,
                'field': METADATA[self.run]['field'],
                'ra': ra,
                'dec': dec,
                'ra1': ra1,
                'ra2': ra2,
                'dec1': dec1,
                'dec2': dec2,
                'zeropoint': self.zeropoint,
                'exptime': self.exptime,
                'time': str(self.hdu.header['DATE-OBS']).encode('ascii'),
                }
        return meta


###########
# FUNCTIONS
###########

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
    metadata = []
    for band in ['halpha', 'r', 'i']:
        log.info('Starting with band {0}'.format(band))
        # Make sure the output directory exists
        util.setup_dir(os.path.join(IMAGE_DESTINATION, band))
        # Retrieve the list of runs
        if band == 'halpha':
            idx_band = 'ha'
        else:
            idx_band = band
        # [constants.IPHASQC_COND_RELEASE]
        runs = constants.IPHASQC['run_'+idx_band]
        # Prepare each run
        result = clusterview.map(prepare_one, runs, block=True)
        metadata.extend(result)

    # Write the metadata to a table
    mycolumns = (str('filename'), str('band'), str('dr2'), str('run'),
                 str('ccd'), str('field'), str('ra'), str('dec'), 
                 str('ra1'), str('ra2'), str('dec1'), str('dec2'),
                 str('zeropoint'), str('exptime'),
                 str('time'))
    rows = list(itertools.chain.from_iterable(metadata)) # flatten list
    t = table.Table(rows, names=mycolumns)
    table_filename=os.path.join(IMAGE_DESTINATION, 'iphas-images.fits')
    t.write(table_filename, format='fits', overwrite=True)





##############################
# MAIN EXECUTION FOR DEBUGGING
##############################

if __name__ == '__main__':
    """
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
    """
    prepare_one(571408)
