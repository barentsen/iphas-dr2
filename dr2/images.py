#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Prepare IPHAS images for data release.

This script will edit all images in the data release to ensure they have
the correct astrometric solution (WCS) and calibration information (zeropoint)
in the header.

TODO
----
- Compress with gzip rather than rice, for wide adoption?

"""
from __future__ import division, print_function, unicode_literals
from astropy.io import fits
from astropy.io import ascii
from astropy import log
import numpy as np
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
MYDESTINATION = os.path.join(constants.DESTINATION, 'images')
util.setup_dir(MYDESTINATION)

# Table containing slight updates to WCS astrometric parameters
WCSFIXES_PATH = os.path.join(constants.PACKAGEDIR, 'wcs-tuning', 'wcs-fixes.csv')
WCSFIXES = ascii.read(WCSFIXES_PATH)

MD = fits.getdata(os.path.join(constants.DESTINATION, 'metadata.fits'))
METADATA = dict(zip(MD['run'], MD))


###########
# CLASSES
###########

"""
class RunDatabase(object):
    def __init__(self):
        table = ascii.read(os.path.join(constants.DESTINATION, 'runs.csv'))
        self.runs = dict(zip(table['run'], table))
"""

class CalibrationDatabase(object):
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

    def __init__(self, run, ccd):
        self.run = run
        self.ccd = ccd
        # Open the image
        self.path_orig = constants.RAWDATADIR + METADATA[run]['image']
        self.fits_orig = fits.open(self.path_orig, do_not_scale_image_data=True)
        # Sort out the new FITS image and header
        self.hdu = fits.PrimaryHDU(self.fits_orig[self.ccd].data)
        self.configure_header()

    @property
    def output_filename(self):
        """Filename of the output?"""
        return '{0}_{1}.fits'.format(self.run, self.ccd)

    @property
    def exptime(self):
        return METADATA[self.run]['exptime_precalib']

    @property
    def zeropoint(self):
        try:
            shift = CALDB.shifts[self.run]
        except KeyError:
            shift = 0.0
        return METADATA[self.run]['zeropoint_precalib'] - self.percorr + shift

    @property
    def percorr(self):
        return METADATA[self.run]['CCD{0}_PERCORR'.format(self.ccd)]

    def configure_header(self):
        # Copy keywords from the original HDU[0]
        for kw in ['RUN', 'OBSERVAT', 'OBSERVER', 'OBJECT',
                   'LATITUDE', 'LONGITUD', 'HEIGHT', 'SLATEL',
                   'TELESCOP', 'RA', 'DEC',
                   'MJD-OBS', 'JD', 'AZIMUTH', 'ZD', 'PLATESCA', 'TELFOCUS',
                   'ROTTRACK', 'ROTSKYPA', 'PARANGLE', 'DOMEAZ', 'AIRMASS',
                   'TEMPTUBE', 'INSTRUME', 'WFFPOS', 'WFFBAND', 'WFFID',
                   'SECPPIX', 'UTSTART', 'DATE-OBS', 'DETECTOR', 'CCDSPEED',
                   'CCDXBIN', 'CCDYBIN', 'CCDSUM', 'CCDTEMP', 'NWINDOWS']:
            self.hdu.header[kw] = self.fits_orig[0].header[kw]
            self.hdu.header.comments[kw] = self.fits_orig[0].header.comments[kw]

        # Copy keywords from the original image extension
        for kw in ['BSCALE', 'BZERO', 'CCDNAME', 'CCDXPIXE', 'CCDYPIXE', 
                   'AMPNAME', 'GAIN', 'READNOIS', 'ZEROCOR', 'LINCOR', 
                   'FLATCOR', 'SEEING', 'WCSPASS', 'NUMBRMS', 'STDCRMS',
                   'PERCORR', 'EXTINCT']:
            self.hdu.header[kw] = self.fits_orig[self.ccd].header[kw]
            self.hdu.header.comments[kw] = self.fits_orig[self.ccd].header.comments[kw]

        # Fix WCS - see fix_wcs() in detections.py!
        # Enforce the pipeline's defaults
        self.hdu.header['RADESYSa'] = 'ICRS'
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

        # Is an updated (fixed) WCS available?
        if self.run in WCSFIXES['RUN']:
            for ccd in constants.EXTENSIONS:
                idx = ((WCSFIXES['RUN'] == self.run)
                       & (WCSFIXES['CCD'] == self.ccd))
                if idx.sum() > 0:
                    log.info("WCS fixed: {0}[{1}].".format(self.run, self.ccd))
                    idx_fix = idx.nonzero()[0][-1]
                    for kw in ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                               'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
                        self.hdu.header[kw] = WCSFIXES[kw][idx_fix]

        
        # Fix zeropoint
        self.hdu.header['ORIGZPT'] = self.fits_orig[self.ccd].header['MAGZPT']
        self.hdu.header.comments['ORIGZPT'] = 'Original nightly zeropoint'
        
        self.hdu.header['MAGZPT'] = self.zeropoint
        self.hdu.header.comments['MAGZPT'] = 'Re-calibrated DR2 zeropoint'

        # Fix exposure time -- it might have changed in detections.py
        self.hdu.header['EXPTIME'] = self.exptime
        self.hdu.header.comments['EXPTIME'] = '[sec] Exposure time assumed by DR2 pipeline'

        # Add history
        for line in str(self.fits_orig[self.ccd].header['HISTORY']).split('\n'):
            self.hdu.header['HISTORY'] = line
        self.hdu.header['HISTORY'] = datetime.datetime.now().strftime('%Y%m%d %H:%M:%S')
        self.hdu.header['HISTORY'] = '    Headers updated by Geert Barentsen as part of DR2.'
        self.hdu.header['HISTORY'] = '    This included changes to MAGZPT, EXPTIME and the WCS.'
        
        self.hdu.header['COMMENT'] = 'Calibration info'
        self.hdu.header['COMMENT'] = '================'
        self.hdu.header['COMMENT'] = 'The MAGZPT keyword in this header has been corrected for atmospheric'
        self.hdu.header['COMMENT'] = 'extinction and gain (PERCORR) and has been re-calibrated as part of DR2.'
        self.hdu.header['COMMENT'] = ''
        self.hdu.header['COMMENT'] = 'Hence to obtain calibrated magnitudes relative to Vega, use:'
        self.hdu.header['COMMENT'] = '    mag(Vega) = MAGZPT - 2.5*log(pixel value / EXPTIME)'

    def save(self, directory=MYDESTINATION):
        target = os.path.join(directory, self.output_filename)
        #prihdu = fits.PrimaryHDU()
        #hdu = fits.HDUList( [prihdu, self.fits[self.ccd]] )
        #hdu.writeto(target, clobber=True)
        self.hdu.writeto(target, clobber=True)

    def summary(self):
        d = {'filename': self.output_filename,
             'exptime': self.exptime,
             'zeropoint': self.zeropoint,
             }
        return d


###########
# FUNCTIONS
###########

def verify_images():
    """Run through all images produced to verify headers."""
    pass


def prepare_one(run):
    result = []
    for ccd in constants.EXTENSIONS:
        img = SurveyImage(run, ccd)
        img.save()
        result.append(img.summary())
    return result

def prepare_images():
    print(prepare_one(358844))


def write_image_table():
    """Produces a table containing the information of images (ccd per ccd).

    filename,ra,dec,zeropoint,percor,exptime,ra1,dec1,ra2,dec2
    """
    pass


##############################
# MAIN EXECUTION FOR DEBUGGING
##############################

if __name__ == '__main__':
    #write_postcalibration_zeropoints()
    prepare_images()
