# -*- coding: utf-8 -*-
"""
This script reads in catalogues produces by the 'imcore' image detection tool
of the Cambridge Astronomical Survey Unit (CASU) pipeline, and transforms
these catalogues into a more user-friendly format which is known in the
data release as the 'iphasDetection' table.

The transformation involves converting fluxes to magnitudes, computing
celestial coordinates and adding flags to indicate quality problems.

Executing this script will also produce a table called 'iphasRuns.csv' which
serves as an index to the available exposures.

Author: Geert Barentsen

TODO
- Take confidence map into account
"""

from astropy.io import fits
from astropy import wcs
import numpy as np
import logging
import os
import datetime
from multiprocessing import Pool


################################
# CONSTANTS & CONFIGURATION
################################

if os.uname()[1] == 'uhppc11.herts.ac.uk':
    # Where are the pipeline-reduced catalogues?
    DATADIR = "/media/0133d764-0bfe-4007-a9cc-a7b1f61c4d1d/iphas"
    # Where to write the output catalogues?
    DESTINATION = "/home/gb/tmp/iphas-dr2/iphasDetection"
else:
    DATADIR = "/car-data/gb/iphas"
    DESTINATION = "/car-data/gb/iphas-dr2/iphasDetection"

# Yale Bright Star Catalogue (Vizier V50), filtered for IPHAS area and V < 4.5
BRIGHTCAT = fits.getdata('lib/BrightStarCat-iphas.fits', 1)

EXTS = [1, 2, 3, 4]  # Which extensions to expect in the fits catalogues?
PXSCALE = 0.333  # Arcsec/pix of the CCD

# Which are the possible filenames of the confidence maps?
CONF_NAMES = {'Halpha': ['Ha_conf.fits', 'Ha_conf.fit',
                         'Halpha_conf.fit',
                         'ha_conf.fits', 'ha_conf.fit',
                         'h_conf.fits', 'h_conf.fit',
                         'Halpha:197_iphas_aug2003_cpm.fit',
                         'Halpha:197_iphas_sep2003_cpm.fit',
                         'Halpha:197_iphas_oct2003_cpm.fit',
                         'Halpha:197_iphas_nov2003_cpm.fit',
                         'Halpha:197_nov2003b_cpm.fit',
                         'Halpha:197_dec2003_cpm.fit',
                         'Halpha:197_jun2004_cpm.fit',
                         'Halpha:197_iphas_jul2004a_cpm.fit',
                         'Halpha:197_iphas_jul2004_cpm.fit',
                         'Halpha:197_iphas_aug2004a_cpm.fit',
                         'Halpha:197_iphas_aug2004b_cpm.fit',
                         'Halpha:197_iphas_dec2004b_cpm.fit'],
              'r': ['r_conf.fit', 'r_conf.fits',
                    'r:214_iphas_aug2003_cpm.fit',
                    'r:214_dec2003_cpm.fit',
                    'r:214_iphas_nov2003_cpm.fit',
                    'r:214_nov2003b_cpm.fit',
                    'r:214_iphas_sep2003_cpm.fit',
                    'r:214_iphas_aug2004a_cpm.fit',
                    'r:214_iphas_aug2004b_cpm.fit',
                    'r:214_iphas_jul2004a_cpm.fit',
                    'r:214_iphas_jul2004_cpm.fit',
                    'r:214_jun2004_cpm.fit'],
              'i': ['i_conf.fit', 'i_conf.fits',
                    'i:215_iphas_aug2003_cpm.fit',
                    'i:215_dec2003_cpm.fit',
                    'i:215_iphas_nov2003_cpm.fit',
                    'i:215_nov2003b_cpm.fit',
                    'i:215_iphas_sep2003_cpm.fit',
                    'i:215_iphas_aug2004a_cpm.fit',
                    'i:215_iphas_aug2004b_cpm.fit',
                    'i:215_iphas_jul2004a_cpm.fit',
                    'i:215_iphas_jul2004_cpm.fit',
                    'i:215_jun2004_cpm.fit']}

# Cache dict to hold the confidence maps for each filter/directory
confmaps = {'Halpha': {}, 'r': {}, 'i': {}}

# Ignore log of negative fluxes
np.seterr(invalid='ignore', divide='ignore')


###############
# CLASSES
###############


class CatalogueException(Exception):
    """ 
    Exception raised when a catalogue has a *known* problem which cannot be
    fixed.
    """
    pass


class DetectionCatalogue():
    """
    Reads in a detection catalogue in the format produced by the Cambridge
    Astronomical Survey Unit (CASU) and transforms it into a UKIDSS-style
    catalogues with user-friendly coordinates, magnitudes and flags.
    """

    def __init__(self, path):
        """ Constructor """
        self.path = path
        self.directory = '/'.join(path.split('/')[:-1])
        self.filename = path.split('/')[-1]
        try:
            self.fits = fits.open(self.path)
        except IOError, e:
            raise CatalogueException('IOError: %s' % e)
        self.validate()  # Raises CatalogueException for dodgy catalogues
        self.cat_path = self.strip_basedir(path)
        self.image_path = self.get_image_path()
        self.conf_path = self.get_conf_path()

    def hdr(self, field, ext=1):
        """Return the value of the header keyword from extension `ext`."""
        return self.fits[ext].header.get(field)

    def validate(self):
        """
        Test the catalogue for known problems, and fix if possible.

        If the catalogue is not suitable for the IPHAS data release,
        a CatalogueException is raised.
        """
        # The OBJECT keyword must start with the word "intphas" or "iphas"
        if not (self.hdr('OBJECT').startswith('intphas')
                or self.hdr('OBJECT').startswith('iphas')):
            raise CatalogueException('Not an IPHAS run, OBJECT = %s' %
                                     self.hdr('OBJECT'))

        # The filter must be one of Halpha/r/i
        if not self.hdr('WFFBAND') in ['Halpha', 'r', 'i']:
            raise CatalogueException('Unexpected filter, WFFBAND = %s' %
                                     self.hdr('WFFBAND'))

        for ccd in EXTS:

            # Early versions of CASU catalogues chave multiple columns 'Blank'
            # Numpy will throw an exception if multiple columns have the same
            # name, so we need to rename these columns.
            n_columns = len(self.fits[ccd].columns)
            for col in range(26, n_columns, 1):
                name = self.fits[ccd].columns[col].name
                if name == 'Blank':
                    self.fits[ccd].columns[col].name = 'Blank%d' % col

            # In early catalogues, the "Number" (SeqNo) field is called "No."
            if self.fits[ccd].columns[0].name == 'No.':
                self.fits[ccd].columns[0].name = 'Number'

            # The headers contain a combination of old- and modern-
            # style WCS parameters for the ZPN projection coefficients, which
            # confuses libwcs. Moreover, in a few cases the keyword values
            # are plainly wrong. Hence we remove the keywords.
            for kw in ['PV1_0', 'PV1_1', 'PV1_2', 'PV1_3',
                       'PV2_0', 'PV2_1', 'PV2_2', 'PV2_3',
                       'PV3_0', 'PV3_1', 'PV3_3', 'PV3_3',
                       'PROJP1', 'PROJP3', 'WAT1_001', 'WAT2_001']:
                del self.fits[ccd].header[kw]

            # ..and enforce the parameters wich have been used by the pipeline
            self.fits[ccd].header['EQUINOX'] = 2000.0
            self.fits[ccd].header['PV2_1'] = 1.0
            self.fits[ccd].header['PV2_3'] = 220.0
            self.fits[ccd].header['CUNIT1'] = 'deg'
            self.fits[ccd].header['CUNIT2'] = 'deg'

            # Some runs do not have date/time stored due to a glitch in the
            # Telescope Control System. We consider this a show-stopper.
            if not 'UTSTART' in self.fits[ccd].header:
                raise CatalogueException('UTSTART keyword missing')
            if not 'DATE-OBS' in self.fits[ccd].header:
                raise CatalogueException('DATE-OBS keyword missing')
            if not 'MJD-OBS' in self.fits[ccd].header:
                raise CatalogueException('MJD-OBS keyword missing')

    def get_image_path(self):
        """Returns the filename of the accompanying image FITS file.

        Raises a CatalogueException if the image is missing.
        """
        candidate = os.path.join(self.directory,
                                 self.filename.split('_')[0] + '.fit')
        if os.path.exists(candidate):
            return self.strip_basedir(candidate)
        else:
            raise CatalogueException('No image found for %s' % (self.path))

    def get_conf_path(self):
        """Return the filename of the accompanying confidence map."""
        mydir = self.directory
        myband = self.hdr('WFFBAND')
        global confmaps
        # The result from previous function calls are stored in 'confmaps'
        if mydir not in confmaps[myband].keys():
            # Some directories do not contain confidence maps
            if mydir == DATADIR+'iphas_nov2006c':
                candidatedir = DATADIR+'iphas_nov2006b'
            elif mydir == DATADIR+'iphas_jul2008':
                candidatedir = DATADIR+'iphas_aug2008'
            elif mydir == DATADIR+'iphas_oct2009':
                candidatedir = DATADIR+'iphas_nov2009'
            elif mydir == DATADIR+'run10':
                candidatedir = DATADIR+'run11'
            elif mydir == DATADIR+'run13':
                candidatedir = DATADIR+'run12'
            else:
                candidatedir = mydir

            # Try all possible names
            for name in CONF_NAMES[myband]:
                candidate = os.path.join(candidatedir, name)
                if os.path.exists(candidate):
                    confmaps[myband][mydir] = candidate  # Success!
                    continue

        # Return confidence map name if we found one, raise exception otherwise
        try:
            return self.strip_basedir(confmaps[myband][mydir])
        except KeyError:
            raise CatalogueException('No confidence map found in %s' % mydir)

    def strip_basedir(self, path):
        return path[len(DATADIR):]

    def get_exptime(self):
        """Return the exposure time for the catalogue, while taking into
        account the quirky nature of EXPTIME as recorded by WFC.

        This follows the original perl script from Brent:

        {
        my ($time) = @_;
           if($time < 15 && abs($time-10) > 0.1){
              return 10.00;
           } elsif ($time > 15 && $time < 35 && abs($time-30) > 0.1){
              return 30.00;
           } elsif ($time > 100 && abs($time-120) > 0.1){
              return 120.00;
           } else {
              return $time;
           }
        }
        """
        t = self.hdr('EXPTIME')
        if t < 15 and abs(t-10) > 0.1:
            return 10.00
        elif t > 15 and t < 35 and abs(t-30) > 0.1:
            return 30.00
        elif t > 100 and abs(t-120) > 0.1:
            return 120.00
        else:
            return t

    def compute_magnitudes(self, n_pixels, flux_field, apcor_field):
        """Convert the flux counts to magnitudes.

        According to the header parameters in the specified extension.
        Be aware that APCOR and PERCORR differ on a CCD-by-CCD basis.
        """
        magnitudes = np.array([])
        exptime = self.get_exptime()

        # Loop over the four CCD extensions and compute the magnitudes
        # Note: APCOR differs on a CCD-by-CCD basis!
        for ccd in EXTS:
            # Mag = ZP - 2.5*log10(flux/exptime) - apcor - percorr
            # See http://apm3.ast.cam.ac.uk/~mike/iphas/README.catalogues
            mypercorr = self.hdr('PERCORR', ccd)
            if mypercorr is None:  # PERCORR keyword is sometimes missing
                mypercorr = 0.0

            flux = self.fits[ccd].data.field(flux_field)
            mag = (self.hdr('MAGZPT', ccd)
                   - 2.5 * np.log10(flux / exptime)
                   - self.hdr(apcor_field, ccd)
                   - mypercorr)
            magnitudes = np.concatenate((magnitudes, mag))

        return magnitudes

    def compute_magnitude_errors(self, n_pixels, flux_field, apcor_field):
        """Convert the flux errors to magnitude errors."""
        errors = np.array([])
        for ccd in EXTS:
            # See http://apm3.ast.cam.ac.uk/~mike/iphas/README.catalogues
            flux = self.fits[ccd].data.field(flux_field)
            err_flux = np.sqrt((flux / self.hdr('GAIN', ccd))
                               + n_pixels * (self.hdr('SKYNOISE', ccd)**2.))
            err_mag = (2.5 / np.log(10)) * err_flux / flux
            errors = np.concatenate((errors, err_mag))
        return errors

    def get_magnitude_column(self, name='aperMag2'):
        mynames = {'peakMag': 'peak', 'peakMagErr': 'peak',
                   'aperMag1': 'core1', 'aperMag1Err': 'core1',
                   'aperMag2': 'core', 'aperMag2Err': 'core',
                   'aperMag3': 'core2', 'aperMag3Err': 'core2',
                   'aperMag4': 'core3', 'aperMag4Err': 'core3',
                   'aperMag5': 'core4', 'aperMag5Err': 'core4'}
        aperture = mynames[name]

        if aperture == 'peak':
            # Peak pixel
            n_pixels = 1
            flux_field = 'Peak_height'
            apcor_field = 'APCORPK'
        elif aperture == 'core1':
            # Radius = 1/2 x rcore
            # Corresponds to Apermag1 in mercats
            n_pixels = np.pi * (0.5*self.hdr('RCORE'))**2
            flux_field = 'Core1_flux'
            apcor_field = 'APCOR1'
        elif aperture == 'core':
            # Radius = rcore
            # Corresponds to Apermag2 in mercats
            n_pixels = np.pi * self.hdr('RCORE')**2
            flux_field = 'Core_flux'
            apcor_field = 'APCOR'
        elif aperture == 'core2':
            # Radius = sqrt(2) x rcore
            # Corresponds to Apermag3 in mercats
            n_pixels = np.pi * (np.sqrt(2.0)*self.hdr('RCORE'))**2
            flux_field = 'Core2_flux'
            apcor_field = 'APCOR2'
        elif aperture == 'core3':
            # Radius = 2 x rcore
            # Corresponds to Apermag4 in mercats
            n_pixels = np.pi * (2.0*self.hdr('RCORE'))**2
            flux_field = 'Core3_flux'
            apcor_field = 'APCOR3'
        elif aperture == 'core4':
            # Radius = 2 sqrt(2) x rcore
            # Corresponds to Apermag5 in mercats
            n_pixels = np.pi * (2.0*np.sqrt(2.0)*self.hdr('RCORE'))**2
            flux_field = 'Core4_flux'
            apcor_field = 'APCOR4'
        else:
            raise CatalogueException('Did not understand requested aperture')

        if name.endswith('Err'):
            errors = self.compute_magnitude_errors(n_pixels, flux_field, 
                                                   apcor_field)
            return fits.Column(name=name, format='E',
                               unit='Sigma', array=errors)
        else:
            mag = self.compute_magnitudes(n_pixels, flux_field, apcor_field)
            return fits.Column(name=name, format='E',
                               unit='Magnitude', array=mag)

    def compute_coordinates(self):
        """Returns RA/DEC using the pixel coordinates and the header WCS"""
        ra = np.array([])
        dec = np.array([])
        for ccd in EXTS:
            mywcs = wcs.WCS(self.fits[ccd].header, relax=True)
            myra, mydec = mywcs.wcs_pix2world(
                            self.fits[ccd].data.field('X_coordinate'),
                            self.fits[ccd].data.field('Y_coordinate'),
                            1)

            ra = np.concatenate((ra, myra))
            dec = np.concatenate((dec, mydec))
        return (ra, dec)

    def compute_brightNeighb(self, ra, dec):
        """ Returns an array of boolean flags indicating whether the stars
        are within 10 arcmin of a star brighter than V < 4.5 """
        flags = np.zeros(len(ra), dtype=bool)  # Initialize result array
        # Try all stars in the truncated bright star catalogue (BSC, Yale)
        # which are nearby-ish
        nearbyish = np.abs(dec[0] - BRIGHTCAT.field('_DEJ2000')) < 2.
        for i in np.where(nearbyish)[0]:
            d_ra = ra - BRIGHTCAT.field('_RAJ2000')[i]
            d_dec = dec - BRIGHTCAT.field('_DEJ2000')[i]
            # Approx angular separation (Astronomical Algorithms Eq. 16.2)
            d = np.sqrt((d_ra*np.cos(np.radians(dec)))**2 + d_dec**2)
            # Flag bright neighbours if within 10 arcmin
            flags[d < 10/60.] = True
        return flags

    def get_csv_summary(self):
        """ Returns a CSV-formatted summary line """
        # Average seeing and ellipticity across CCDs
        avg_seeing = (PXSCALE * np.mean([self.hdr('SEEING', ccd)
                                         for ccd in EXTS]))
        avg_ellipt = np.mean([self.hdr('ELLIPTIC', i) for i in EXTS])

        # When the PERCORR keyword is missing, assume it is zero
        mypercorr = self.hdr('PERCORR')
        if mypercorr is None:
            mypercorr = 0.0

        e_5sig = (self.hdr('MAGZPT')
                  - 2.5 * np.log10(
                              5.0*np.sqrt(
                                  np.sqrt(2.0)*np.pi*self.hdr('RCORE')**2.)
                                  * self.hdr('SKYNOISE') / self.get_exptime())
                  - self.hdr('APCOR2')
                  - mypercorr)

        field = self.hdr('OBJECT').split('_')[1].split(' ')[0]

        return ('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,'
                + '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,'
                + '%s,%s,%s,%s,%s,%s,"%s",%s,%s,%s,%s,%s,%s,%s') % (
                    self.cat_path,
                    self.image_path,
                    self.conf_path,
                    self.hdr('RUN'),
                    self.hdr('OBJECT'),
                    self.hdr('RA'),
                    self.hdr('DEC'),
                    field,
                    avg_seeing,
                    PXSCALE*self.hdr('SEEING', 1),
                    PXSCALE*self.hdr('SEEING', 2),
                    PXSCALE*self.hdr('SEEING', 3),
                    PXSCALE*self.hdr('SEEING', 4),
                    avg_ellipt,
                    self.hdr('ELLIPTIC', 1), self.hdr('ELLIPTIC', 2),
                    self.hdr('ELLIPTIC', 3), self.hdr('ELLIPTIC', 4),
                    e_5sig,
                    self.hdr('AIRMASS'),
                    self.hdr('RCORE'),
                    self.hdr('CROWDED'),
                    self.hdr('SKYLEVEL', 1), self.hdr('SKYLEVEL', 2),
                    self.hdr('SKYLEVEL', 3), self.hdr('SKYLEVEL', 4),
                    self.hdr('SKYNOISE', 1), self.hdr('SKYNOISE', 2),
                    self.hdr('SKYNOISE', 3), self.hdr('SKYNOISE', 4),
                    self.hdr('MAGZPT'),
                    self.hdr('MAGZRR'),
                    self.hdr('PERCORR', 1), self.hdr('PERCORR', 2),
                    self.hdr('PERCORR', 3), self.hdr('PERCORR', 4),
                    self.hdr('GAIN', 1), self.hdr('GAIN', 2),
                    self.hdr('GAIN', 3), self.hdr('GAIN', 4),
                    self.hdr('STDCRMS', 1), self.hdr('STDCRMS', 2),
                    self.hdr('STDCRMS', 3), self.hdr('STDCRMS', 4),
                    self.hdr('CCDSPEED'),
                    self.hdr('OBSERVER'),
                    self.hdr('DAZSTART'),
                    self.hdr('DATE-OBS')+' '+self.hdr('UTSTART'),
                    self.hdr('MJD-OBS'),
                    self.hdr('EXPTIME'),
                    self.hdr('WFFPOS'),
                    self.hdr('WFFBAND'),
                    self.hdr('WFFID')
                )

    def concat(self, name):
        """Returns the concatenated array of a column across all ccds"""
        if name in self.fits[1].columns.names:
            return np.concatenate([self.fits[ccd].data.field(name)
                                   for ccd in EXTS])
        else:
            # If the key does not exist, return an array of None's
            return np.concatenate([[None]*self.fits[ccd].data.size
                                   for ccd in EXTS])

    def save_detections(self):
        """ Create the columns of the output FITS table and save them.

        The fits data types used are:
        L = boolean (1 byte?)
        X = bit
        B = unsigned byte
        I = 16-bit int
        J = 32-bit int
        K = 64-bit int
        A = 1-byte char
        E = single
        D = double
        """
        output_filename = os.path.join(DESTINATION,
                                       '%s_det.fits' % self.hdr('RUN'))

        # Total number of detected objects
        n_objects = np.sum([self.fits[ccd].data.size for ccd in EXTS])

        # Pre-prepare columns
        ccds = np.concatenate([[ccd] * self.fits[ccd].data.size
                              for ccd in EXTS])
        runID = np.array([self.hdr('RUN')] * n_objects)
        seqNo = self.concat('Number')
        bandnames = {'r': 'r', 'i': 'i', 'Halpha': 'ha'}
        myband = bandnames[self.hdr('WFFBAND')]
        band = np.array([myband] * n_objects)

        detectionID = np.array([int('%07d%d%06d' % (
                    self.hdr('RUN'), ccds[i], seqNo[i]))
                    for i in range(n_objects)])
        col_detectionID = fits.Column(name='detectionID', format='K',
                                      unit='Number', array=detectionID)
        col_runID = fits.Column(name='runID', format='J',
                                unit='Number', array=runID)
        col_ccd = fits.Column(name='ccd', format='B', unit='Number',
                              array=ccds)
        col_seqNum = fits.Column(name='seqNum', format='I', unit='Number',
                                 array=seqNo)
        col_band = fits.Column(name='band', format='2A', unit='String',
                               array=band)
        col_x = fits.Column(name='x', format='E', unit='Pixels',
                            array=self.concat('X_coordinate'))
        col_y = fits.Column(name='y', format='E', unit='Pixels',
                            array=self.concat('Y_coordinate'))
        # Double precision equatorial coordinates
        myra, mydec = self.compute_coordinates()

        col_ra = fits.Column(name='ra', format='D', unit='deg',
                             array=myra)  # Double precision!
        col_dec = fits.Column(name='dec', format='D', unit='deg',
                              array=mydec)  # Double precision!
        # Astrometric fit error (arcsec)
        posErr = np.concatenate([
                    [self.fits[ccd].header.get('STDCRMS')]
                    * self.fits[ccd].data.size
                    for ccd in EXTS])  # In arcsec
        col_posErr = fits.Column(name='posErr', format='E', unit='arcsec',
                                 array=posErr)  # Double precision!
        col_gauSig = fits.Column(name='gauSig', format='E', unit='Number',
                                 array=self.concat('Gaussian_sigma'))
        col_ell = fits.Column(name='ell', format='E', unit='Number',
                              array=self.concat('Ellipticity'))
        col_pa = fits.Column(name='pa', format='E', unit='Number',
                             array=self.concat('Position_angle'))
        col_sky = fits.Column(name='sky', format='E', unit='Counts',
                              array=self.concat('Skylev'))
        col_skyVar = fits.Column(name='skyVar', format='E', unit='Counts',
                                 array=self.concat('Skyrms'))
        col_class = fits.Column(name='class', format='I', unit='Flag',
                                array=self.concat('Classification'))
        col_classStat = fits.Column(name='classStat', format='E',
                                    unit='N-sigma',
                                    array=self.concat('Statistic'))

        # For deblended images, only the 1st areal profile is computed
        # and the rest are set to -1
        deblend = (self.concat('Areal_3_profile') < 0)
        col_deblend = fits.Column(name='deblend', format='L', unit='Boolean',
                                  array=deblend)

        # Which stars are saturated?
        # The saturation level is stored in the SATURATE keyword for each ccd
        saturated = np.concatenate([(
                        self.fits[ccd].data.field('Peak_height')
                        > self.fits[ccd].header.get('SATURATE'))
                        for ccd in EXTS])
        col_saturated = fits.Column(name='saturated', format='L',
                                    unit='Boolean', array=saturated)

        # Mark stars near the edges
        avoidance = 4.0/0.333  # 4 Arcseconds
        min_x = 1 + avoidance
        max_x = 2048 - avoidance
        min_y = 1 + avoidance
        max_y = 4096 - avoidance

        truncated = ((self.concat('X_coordinate') < min_x)
                     | (self.concat('X_coordinate') > max_x)
                     | (self.concat('Y_coordinate') < min_y)
                     | (self.concat('Y_coordinate') > max_y))
        col_truncated = fits.Column(name='truncated', format='L',
                                    unit='Boolean', array=truncated)

        brightNeighb = self.compute_brightNeighb(myra, mydec)
        col_brightNeighb = fits.Column(name='brightNeighb', format='L',
                                       unit='Boolean', array=brightNeighb)

        badPix = self.concat('Bad_pixels')
        col_badPix = fits.Column(name='badPix', format='E', unit='Pixels',
                                 array=badPix)

        # Figure out the YYYYMMDD identifier of the *night* (i.e. evening)
        mydate = datetime.datetime.strptime(
                        self.hdr('DATE-OBS')+' '+self.hdr('UTSTART')[0:2],
                        '%Y-%m-%d %H')  # Dont parse seconds; they can be '60'
        if mydate.hour < 12:
            mydate -= datetime.timedelta(1)
        night = np.array([mydate.strftime('%Y%m%d')] * n_objects)
        col_night = fits.Column(name='night', format='J',
                                array=night)

        mjd = np.array([self.hdr('MJD-OBS')] * n_objects)
        col_mjd = fits.Column(name='mjd', format='D', unit='Julian days',
                              array=mjd)

        seeing = np.concatenate([[PXSCALE*self.hdr('SEEING', ccd)]
                                 * self.fits[ccd].data.size
                                 for ccd in EXTS])
        col_seeing = fits.Column(name='seeing', format='E', unit='arcsec',
                                 array=seeing)

        # Write the output fits table
        cols = fits.ColDefs([col_detectionID, col_runID,
                             col_ccd, col_seqNum, col_band,
                             col_x, col_y,
                             col_ra, col_dec, col_posErr,
                             col_gauSig, col_ell, col_pa,
                             self.get_magnitude_column('peakMag'),
                             self.get_magnitude_column('peakMagErr'),
                             self.get_magnitude_column('aperMag1'),
                             self.get_magnitude_column('aperMag1Err'),
                             self.get_magnitude_column('aperMag2'),
                             self.get_magnitude_column('aperMag2Err'),
                             self.get_magnitude_column('aperMag3'),
                             self.get_magnitude_column('aperMag3Err'),
                             col_sky, col_skyVar,
                             col_class, col_classStat, col_badPix,
                             col_deblend, col_saturated,
                             col_truncated, col_brightNeighb,
                             col_night, col_mjd, col_seeing])
        hdu_table = fits.new_table(cols, tbtype='BinTableHDU')
        # Copy some of the original keywords
        for kw in ['RUN', 'OBSERVAT', 'LATITUDE', 'LONGITUD', 'HEIGHT',
                   'OBSERVER',
                   'OBJECT', 'RA', 'DEC', 'EQUINOX', 'RADECSYS',
                   'MJD-OBS', 'JD', 'DATE-OBS', 'UTSTART',
                   'INSTRUME', 'WFFPOS', 'WFFBAND', 'WFFPSYS', 'WFFID',
                   'EXPTIME', 'AIRMASS', 'MAGZPT', 'MAGZRR']:
            hdu_table.header[kw] = self.hdr(kw, 4)

        for ext in EXTS:
            hdu_table.header['SEEING%d' % ext] = self.hdr('SEEING', ext)
            hdu_table.header['ELLIP%d' % ext] = self.hdr('ELLIPTIC', ext)
            hdu_table.header['SKY%d' % ext] = self.hdr('SKYLEVEL', ext)
        hdu_table.header['EXPTUSED'] = self.get_exptime()
        hdu_table.header['CATALOG'] = self.cat_path
        hdu_table.header['IMAGE'] = self.image_path
        hdu_table.header['CONFMAP'] = self.conf_path

        hdu_primary = fits.PrimaryHDU()
        hdulist = fits.HDUList([hdu_primary, hdu_table])
        hdulist.writeto(output_filename, clobber=True)


######################
# FUNCTIONS
######################


def list_catalogues():
    logging.info('Searching for catalogues in %s' % DATADIR)
    catalogues = []
    for mydir in os.walk(DATADIR, followlinks=True):
        logging.info('Entering %s' % mydir[0])
        for filename in mydir[2]:
            # Only consider files of the form *_cat.fits
            if filename.endswith("_cat.fits"):
                catalogues.append(os.path.join(mydir[0], filename))
    logging.info('Found %d catalogues' % len(catalogues))
    return catalogues


def run_one(path):
    try:
        pc = DetectionCatalogue(path)
        csv = pc.get_csv_summary()
        pc.save_detections()
        return csv
    except CatalogueException, e:
        logging.error('%s: CatalogueException: %s' % (path, e))
        return None
    except Exception, e:
        logging.error('%s: *UNEXPECTED EXCEPTION*: %s' % (path, e))
        return None


def run_all(ncores=4):
    """
    Create catalogues for all runs found in the data directory.

    ncores: number of processing cores to use
    """
    catalogues = list_catalogues()

    # Execute our analysis for each mercat
    p = Pool(processes=ncores)
    results = p.imap(run_one, catalogues)  # returns an iterator

    # Write the results
    filename = os.path.join(DESTINATION, 'iphasRun.csv')
    out = open(filename, 'w')

    out.write('catalogue,image,conf,run,object,ra,dec,field,'
              + 'seeing,seeing1,seeing2,seeing3,seeing4,'
              + 'ellipt,ellipt1,ellipt2,ellipt3,ellipt4,'
              + '5sig,'
              + 'airmass,rcore,crowded,'
              + 'skylevel1,skylevel2,skylevel3,skylevel4,'
              + 'skynoise1,skynoise2,skynoise3,skynoise4,magzpt,magzrr,'
              + 'percorr1,percorr2,percorr3,percorr4,'
              + 'gain1,gain2,gain3,gain4,'
              + 'stdcrms1,stdcrms2,stdcrms3,stdcrms4,'
              + 'ccdspeed,observer,'
              + 'dazstart,time,mjd-obs,exptime,wffpos,wffband,wffid\n')
    for r in results:
        if r is None:
            continue
        out.write(r+'\n')
        out.flush()
    out.close()


###################
# MAIN EXECUTION
###################

if __name__ == '__main__':
    logging.basicConfig(level=logging.ERROR)
    #run_all(4)

    #Testcases:
    run_one('/media/0133d764-0bfe-4007-a9cc-a7b1f61c4d1d/iphas/iphas_nov2003b/r375399_cat.fits')
    run_one('/media/0133d764-0bfe-4007-a9cc-a7b1f61c4d1d/iphas/iphas_nov2003b/r375400_cat.fits')
    run_one('/media/0133d764-0bfe-4007-a9cc-a7b1f61c4d1d/iphas/iphas_nov2003b/r375401_cat.fits')

    run_one('/media/0133d764-0bfe-4007-a9cc-a7b1f61c4d1d/iphas/iphas_aug2004a/r414199_cat.fits')
