#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Prepares IPHAS images for public release.

This script will edit all images in the data release to ensure they have
the correct astrometric solution (WCS) and calibration information (PHOTZP).
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
        self.path = constants.RAWDATADIR + METADATA[run]['image']
        self.fits = fits.open(self.path, do_not_scale_image_data=True)

        # Is the run a DR2-recalibrated run?
        mycaldb = get_caldb()
        if self.run in mycaldb.shifts:
            self.calibrated = True
        else:
            self.calibrated = False

        # Sort out the new FITS image and header
        self.hdu = self.fits[self.ccd]
        self.set_header()
        self.fix_wcs()
        self.add_comments()

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
                   'LATITUDE', 'LONGITUD', 'HEIGHT', 'SLATEL', 'TELESCOP',
                   'MJD-OBS', 'JD', 'PLATESCA', 'TELFOCUS', 'AIRMASS',
                   'TEMPTUBE', 'INSTRUME', 'WFFPOS', 'WFFBAND', 'WFFID',
                   'SECPPIX', 'DETECTOR', 'CCDSPEED',
                   'CCDXBIN', 'CCDYBIN', 'CCDSUM', 'CCDTEMP', 'NWINDOWS']:
            self.hdu.header.insert('NAXIS2', (kw,
                                              self.fits[0].header[kw],
                                              self.fits[0].header.comments[kw])
                                   )

        # Ensure a proper ISO stamp
        isostamp = (self.fits[0].header['DATE-OBS']
                    + 'T' + self.fits[0].header['UTSTART'])
        self.hdu.header.insert('NAXIS2', ('DATE-OBS',
                                          isostamp,
                                          'Start time of the exposure [UTC]'))

        # Fix exposure time -- it might have changed in detections.py
        self.hdu.header['EXPTIME'] = self.exptime
        self.hdu.header.comments['EXPTIME'] = '[sec] Exposure time adopted in DR2'

        # Add true zeropoint with all corrections absorbed
        self.hdu.header['PHOTZP'] = self.photzp
        self.hdu.header.comments['PHOTZP'] = 'mag(Vega) = -2.5*log(pixel value) + PHOTZP'
        self.hdu.header.comments['MAGZPT'] = 'Uncorrected nightly ZP (per second)'

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

        # Where is the conf map?
        self.hdu.header['CONFMAP'] = self.confmap

    def fix_wcs(self):
        """Derived from fix_wcs() in detections.py."""
        # Never trust these WCS keywords, which may have been left behind
        # by older versions of the CASU pipeline:
        for kw in ['PV1_0', 'PV1_1', 'PV1_2', 'PV1_3',
                   'PV2_0', 'PV2_1', 'PV2_2', 'PV2_3', 'PV2_5',
                   'PV3_0', 'PV3_1', 'PV3_3', 'PV3_3',
                   'PROJP1', 'PROJP3', 'PROJP5', 'WAT1_001', 'WAT2_001',
                   'RADECSYS', 'CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2']:
            try:
                del self.hdu.header[kw]
            except KeyError:
                pass

        # Ensure the pipeline defaults are set correctly
        self.hdu.header.insert('CRVAL1', ('RADESYS', 'ICRS', 'WCS calibrated against 2MASS'))
        self.hdu.header.insert('CRVAL1', ('EQUINOX', 2000.0))
        self.hdu.header.insert('CRVAL1', ('CTYPE1', 'RA---ZPN', 'Algorithm type for axis 1'))
        self.hdu.header.insert('CRVAL1', ('CTYPE2', 'DEC--ZPN', 'Algorithm type for axis 2'))
        self.hdu.header.insert('CRVAL1', ('CRUNIT1', 'deg', 'Unit of right ascension coordinates'))
        self.hdu.header.insert('CRVAL1', ('CRUNIT2', 'deg', 'Unit of declination coordinates'))
        self.hdu.header.insert('CRVAL1', ('PV2_1', 1.0, 'Coefficient for r term'))
        self.hdu.header.insert('CRVAL1', ('PV2_2', 0.0, 'Coefficient for r**2 term'))
        self.hdu.header.insert('CRVAL1', ('PV2_3', 220.0, 'Coefficient for r**3 term'))

        # Improvide the documentation following
        # http://apm49.ast.cam.ac.uk/surveys-projects/wfcam/technical/astrometry
        self.hdu.header.comments['CRPIX1'] = '[pixel] Reference pixel along axis 1'
        self.hdu.header.comments['CRPIX2'] = '[pixel] Reference pixel along axis 2'
        self.hdu.header.comments['CRVAL1'] = '[deg] Right ascension at the reference pixel'
        self.hdu.header.comments['CRVAL2'] = '[deg] Declination at the reference pixel'
        self.hdu.header.comments['CD1_1'] = 'Transformation matrix element'
        self.hdu.header.comments['CD1_2'] = 'Transformation matrix element'
        self.hdu.header.comments['CD2_1'] = 'Transformation matrix element'
        self.hdu.header.comments['CD2_2'] = 'Transformation matrix element'

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
        """Populate the HISTORY and COMMENT keywords of the FITS file.

        At the time of writing, we use "hdu._header" rather than "hdu.header"
        to by-pass Astropy issue #2363.
        """
        self.hdu._header['HISTORY'] = datetime.datetime.now().strftime('%Y-%m-%d')+':'
        self.hdu._header['HISTORY'] = 'Headers updated by G. Barentsen to add the DR2 calibration information.'
        self.hdu._header['HISTORY'] = 'This includes changes to PHOTZP, EXPTIME, and the WCS keywords.'

        # Set calibration comments
        if self.calibrated:
            self.hdu._header['COMMENT'] = 'Photometric calibration info'
            self.hdu._header['COMMENT'] = '============================'
            self.hdu._header['COMMENT'] = 'The pixel values (number counts) in this image can be converted into'
            self.hdu._header['COMMENT'] = 'Vega-based magnitudes using the PHOTZP keyword as follows:'
            self.hdu._header['COMMENT'] = '    mag(Vega) = -2.5*log10(pixel value) + PHOTZP.'
            self.hdu._header['COMMENT'] = 'The PHOTZP value has been computed such that it absorbs the required'
            self.hdu._header['COMMENT'] = 'corrections for atmospheric extinction, gain variations, exposure time,'
            self.hdu._header['COMMENT'] = 'and the DR2 re-calibration shift.'
            self.hdu._header['COMMENT'] = 'As these images still include moonlight and other sources of'
            self.hdu._header['COMMENT'] = 'non-astronomical background, they can only support flux measurements'
            self.hdu._header['COMMENT'] = 'that include a suitably-chosen local background subtraction.'
        else:
            self.hdu._header['COMMENT'] = 'WARNING'
            self.hdu._header['COMMENT'] = '======='
            self.hdu._header['COMMENT'] = 'This image is not part of IPHAS DR2 and has not been re-calibrated.'
            self.hdu._header['COMMENT'] = 'It may have been excluded from DR2 due to a serious quality problem,'
            self.hdu._header['COMMENT'] = 'for example, PHOTZP may be inaccurate due to clouds.'
            self.hdu._header['COMMENT'] = '*USE AT YOUR OWN RISK*'

    @property
    def output_filename(self):
        """Filename of the output?"""
        return 'r{0}-{1}.fits.fz'.format(self.run, self.ccd).encode('ascii')

    def save(self):
        """Save the ccd image to a new file."""
        directory = os.path.join(constants.PATH_IMAGES,
                                 'r'+str(self.run)[0:3])
        util.setup_dir(directory)
        target = os.path.join(directory, self.output_filename)
        # checksum=True adds the CHECKSUM and DATASUM keywords
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
                'in_dr2': in_dr2,
                'ra': ra,
                'dec': dec,
                'ra_min': ra1,
                'ra_max': ra2,
                'dec_min': dec1,
                'dec_max': dec2,
                'band': band,
                'utstart': str(self.hdu.header['DATE-OBS']).encode('ascii'),
                'exptime': self.hdu.header['EXPTIME'],
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
                img.fits.close()  # avoid memory leak
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
        result = clusterview.map(prepare_one, runs, block=True)
        metadata.extend(result)

    # Write the metadata to a table
    mycolumns = (str('filename'), str('run'), str('ccd'),
                 str('in_dr2'),
                 str('ra'), str('dec'),
                 str('ra_min'), str('ra_max'),
                 str('dec_min'), str('dec_max'),
                 str('band'),
                 str('utstart'), str('exptime'),
                 str('seeing'), str('elliptic'),
                 str('skylevel'), str('skynoise'),
                 str('airmass'), str('photzp'),
                 str('confmap'))
    rows = list(itertools.chain.from_iterable(metadata))  # flatten list
    t = table.Table(rows, names=mycolumns)
    table_filename = os.path.join(constants.PATH_IMAGES, 'iphas-images.fits')
    t.write(table_filename, format='fits', overwrite=True)


################################
# MAIN EXECUTION (FOR DEBUGGING)
################################

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
    #prepare_one(571408)
    #prepare_one(476845)
