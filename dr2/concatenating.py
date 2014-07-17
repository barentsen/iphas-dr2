#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Concatenates the seamed and bandmerged fields into the final catalogue.

This script takes the output from the seaming script and concatenates the
results into the final source catalogue products, which contains one entry
for each unique source and is generated in 5x5 degree tiles. Both a 'light' and 
a 'full' version of these tiles are generated.

This step also creates the source designations "JHHMMSS.ss+DDMMSS.s".
"""
from __future__ import division, print_function, unicode_literals
import os
import numpy as np
import datetime
from multiprocessing import Pool
from astropy import log

import constants
from constants import IPHASQC
import util

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Geert Barentsen', 'Hywel Farnhill', 'Janet Drew']


###########
# CLASSES
###########

class Concatenator(object):
    """Concatenates field-based data into a partial catalogue."""

    def __init__(self, strip, part='a', mode='full'):
        assert(part in ['a', 'b'])
        assert(mode in ['light', 'full'])

        self.strip = strip
        self.part = part
        self.mode = mode
        
        # Where are the input catalogues?
        self.datapath = os.path.join(constants.DESTINATION, 'seamed')

        # Where to write the output?
        self.destination = os.path.join(constants.DESTINATION,
                                        'concatenated')
        
        # Setup the destination directory
        if mode == 'light':
            self.destination = os.path.join(self.destination, 'light')
        else:
            self.destination = os.path.join(self.destination, 'full')
        util.setup_dir(self.destination)
        util.setup_dir(self.destination+'-compressed')

        log.info('Reading data from {0}'.format(self.datapath))

        # Limits
        self.lon1 = strip
        self.lon2 = strip + constants.STRIPWIDTH
        self.fieldlist = self.get_fieldlist()

    def get_partname(self):
        """Returns the name of this partial catalogue, e.g. '215b'"""
        return '{0:03.0f}{1}'.format(self.lon1, self.part)

    def get_output_filename(self, gzip=False):
        """Returns the full path of the output file."""
        if self.mode == 'light':
            suffix = '-light'
        else:
            suffix = ''

        destination = self.destination
        extension = 'fits'
        if gzip:
            destination += '-compressed'
            extension += '.gz'
        return os.path.join(destination,
                            'iphas-dr2-{0}{1}.{2}'.format(
                                                    self.get_partname(),
                                                    suffix,
                                                    extension))

    def get_fieldlist(self):
        # Which are our fields?
        # Note: we must allow for border overlaps
        if self.part == 'a':
            cond_b = IPHASQC['b'] < (0 + constants.FIELD_MAXDIST)
        else:
            cond_b = IPHASQC['b'] > (0 - constants.FIELD_MAXDIST)

        cond_strip = (constants.IPHASQC_COND_RELEASE
                      & cond_b
                      & (IPHASQC['l'] >= (self.lon1 - constants.FIELD_MAXDIST))
                      & (IPHASQC['l'] < (self.lon2 + constants.FIELD_MAXDIST)))
        fieldlist = IPHASQC['id'][cond_strip]
        log.info('Found {0} fields.'.format(len(fieldlist)))
        return fieldlist

    def run(self):
        """Performs the concatenation of the strip.

        This step will only keep stars with low errbits: 
            (errBits < 64)
        and not uber-saturated:
            ! (r<12.5 & i<11.5 & ha<12)
        and reasonable errors: 
            (rErr < 0.198 || iErr < 0.198 || haErr < 0.198)
        and not noise-like:
            (pStar > 0.2 || pGalaxy > 0.2)
        and not saturated in all bands:
            (NULL_rErrBits || NULL_iErrBits || NULL_haErrBits || ((rErrbits & iErrBits & haErrBits & 8) == 0))
        and being the primary detection:
            sourceID == primaryID
        """
        if self.part == 'a':
            cond_latitude = "b < 0"
        else:
            cond_latitude = "b >= 0"

        if self.mode == 'full':
            extracmd = """delcols "pSaturated \
                                   rErrBits iErrBits haErrBits errBits \
                                   rPlaneX rPlaneY iPlaneX iPlaneY \
                                   haPlaneX haPlaneY rAxis primaryID \
                                   vignetted truncated badPix" """
        else:
            # select "nBands == 3"; \
            extracmd = """keepcols "name ra dec \
                                    r rErr \
                                    i iErr \
                                    ha haErr \
                                    mergedClass errBits";"""

        instring = ''
        for field in self.fieldlist:
            path = os.path.join(self.datapath,
                                'strip{0:.0f}'.format(self.strip),
                                '{0}.fits'.format(field))
            instring += 'in={0} '.format(path)

        output_filename = self.get_output_filename()
        output_filename_gzip = self.get_output_filename(gzip=True)
        log.info('Writing data to {0}'.format(output_filename))

        version = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')

        # A bug in stilts causes long fieldIDs to be truncated if -utype S15 is not set
        # We also replace a bunch of column descriptions because they cannot be longer than 73 chars.
        param = {'stilts': constants.STILTS,
                 'in': instring,
                 'icmd': """'clearparams *; \
                             setparam NAME "IPHAS DR2 Source Catalogue (part """+self.get_partname()+""")"; \
                             setparam ORIGIN "www.iphas.org"; \
                             setparam AUTHOR "Geert Barentsen, Hywel Farnhill, Janet Drew"; \
                             setparam VERSION \""""+version+""""; \
                             select "(errBits < 64) \
                                      & ! (r<12.5 & i<11.5 & ha<12) \
                                      & (rErr < 0.198 || iErr < 0.198 || haErr < 0.198) \
                                      & (pStar > 0.2 || pGalaxy > 0.2) \
                                      & (NULL_rErrBits || NULL_iErrBits || NULL_haErrBits || ((rErrbits & iErrBits & haErrBits & 8) == 0))
                                      & l >= """+str(self.lon1)+""" \
                                      & l < """+str(self.lon2)+""" \
                                      & """+str(cond_latitude)+""" \
                                      & sourceID == primaryID"; \
                             addcol -before ra \
                                    -desc "Source designation (JHHMMSS.ss+DDMMSS.s) without IPHAS2 prefix." \
                                    name \
                                    "concat(\\"J\\", 
                                            replaceAll(degreesToHms(ra, 2),
                                                       \\":\\", \\"\\"), 
                                            replaceAll(degreesToDms(dec, 1),
                                                       \\":\\", \\"\\")
                                            )"; \
                             addcol -before rMJD -desc "True if source was blended with a nearby neighbour in the r-band." \
                                    rDeblend "NULL_rErrBits ? false : (rErrBits & 2) > 0";
                             addcol -before rMJD -desc "True i the peak pixel count exceeded 55000 in r." \
                                    rSaturated "r<13 ? true : NULL_rErrBits ? false : (rErrBits & 8) > 0";
                             addcol -before iMJD -desc "True if source was blended with a nearby neighbour in the i-band." \
                                    iDeblend "NULL_iErrBits ? false : (iErrBits & 2) > 0";
                             addcol -before iMJD -desc "True if the peak pixel count exceeded 55000 in i." \
                                    iSaturated "i<12 ? true : NULL_iErrBits ? false : (iErrBits & 8) > 0";
                             addcol -before haMJD -desc "True if source was blended with a nearby neighbour in H-alpha." \
                                    haDeblend "NULL_haErrBits ? false : (haErrBits & 2) > 0";
                             addcol -before haMJD -desc "True if the peak pixel count exceeded 55000 in H-alpha." \
                                    haSaturated "ha<12.5 ? true : NULL_haErrBits ? false : (haErrBits & 8) > 0";
                             replacecol saturated "rSaturated || iSaturated || haSaturated";
                             colmeta -name a10 reliable;
                             replacecol a10 "! saturated & nBands == 3 & rErr<0.1 & iErr<0.1 & haErr<0.1 & (abs(r-rAperMag1) < 3*hypot(rErr,rAperMag1Err)+0.03) & (abs(i-iAperMag1) < 3*hypot(iErr,iAperMag1Err)+0.03) & (abs(ha-haAperMag1) < 3*hypot(haErr,haAperMag1Err)+0.03)";
                             addcol -before fieldID -desc "True if (a10 & pStar > 0.9 & ! deblend & ! brightNeighb)" \
                                    a10point "a10 & pStar > 0.9 & ! deblend & ! brightNeighb";
                             replacecol -utype S15 fieldID "fieldID";
                             replacecol -utype S1 fieldGrade "toString(fieldGrade)";
                             colmeta -desc "True if detected in all bands at 10-sigma plus other criteria." a10;
                             colmeta -desc "J2000 RA with respect to the 2MASS reference frame." ra;
                             colmeta -desc "Unique source identification string (run-ccd-detectionnumber)." sourceID;
                             colmeta -desc "Astrometric fit error (RMS) across the CCD." posErr;
                             colmeta -desc "1=galaxy, 0=noise, -1=star, -2=probableStar, -3=probableGalaxy." mergedClass;
                             colmeta -desc "N(0,1) stellarness-of-profile statistic." mergedClassStat;
                             colmeta -desc "1=galaxy, 0=noise, -1=star, -2=probableStar, -3=probableGalaxy." rClass;
                             colmeta -desc "1=galaxy, 0=noise, -1=star, -2=probableStar, -3=probableGalaxy." iClass;
                             colmeta -desc "1=galaxy, 0=noise, -1=star, -2=probableStar, -3=probableGalaxy." haClass;
                             colmeta -desc "Unique r-band detection identifier (run-ccd-detectionnumber)." rDetectionID;
                             colmeta -desc "Unique i-band detection identifier (run-ccd-detectionnumber)." iDetectionID;
                             colmeta -desc "Unique H-alpha detection identifier (run-ccd-detectionnumber)." haDetectionID;
                             colmeta -desc "CCD pixel coordinate in the r-band exposure." rX;
                             colmeta -desc "CCD pixel coordinate in the r-band exposure." rY;
                             colmeta -desc "CCD pixel coordinate in the i-band exposure." iX;
                             colmeta -desc "CCD pixel coordinate in the i-band exposure." iY;
                             colmeta -desc "CCD pixel coordinate in the H-alpha exposure." haX;
                             colmeta -desc "CCD pixel coordinate in the H-alpha exposure." haY;
                             colmeta -desc "Survey field identifier." fieldID;
                             colmeta -desc "Probability the source is extended." pGalaxy;
                             colmeta -desc "Default r mag (Vega) using the 2.3 arcsec aperture." r;
                             colmeta -desc "Default i mag (Vega) using the 2.3 arcsec aperture." i;
                             colmeta -desc "Default H-alpha mag (Vega) using the 2.3 arcsec aperture." ha;
                             colmeta -desc "r mag (Vega) derived from peak pixel height." rPeakMag;
                             colmeta -desc "i mag (Vega) derived from peak pixel height." iPeakMag;
                             colmeta -desc "H-alpha mag (Vega) derived from peak pixel height." haPeakMag;
                             colmeta -desc "r mag (Vega) using the 1.2 arcsec aperture." rAperMag1;
                             colmeta -desc "i mag (Vega) using the 1.2 arcsec aperture." iAperMag1;
                             colmeta -desc "H-alpha mag (Vega) using the 1.2 arcsec aperture." haAperMag1;
                             colmeta -desc "r mag (Vega) using the 3.3 arcsec aperture." rAperMag3;
                             colmeta -desc "i mag (Vega) using the 3.3 arcsec aperture." iAperMag3;
                             colmeta -desc "H-alpha mag (Vega) using the 3.3 arcsec aperture." haAperMag3;
                             colmeta -desc "Internal quality control score of the field. One of A, B, C or D." fieldGrade;
                             colmeta -desc "Number of repeat observations of this source in the survey." nObs;
                             colmeta -desc "SourceID of the object in the partner exposure." sourceID2;
                             colmeta -desc "FieldID of the partner detection." fieldID2;
                             colmeta -desc "r mag (Vega) in the partner field, obtained within 10 minutes." r2;
                             colmeta -desc "Uncertainty for r2." rErr2;
                             colmeta -desc "i mag (Vega) in the partner field, obtained within 10 minutes." i2;
                             colmeta -desc "Uncertainty for i2." iErr2;
                             colmeta -desc "H-alpha mag (Vega) in the partner field, obtained within 10 minutes." ha2;
                             colmeta -desc "Uncertainty for ha2." haErr2;
                             colmeta -desc "flag brightNeighb (1), deblend (2), saturated (8), vignetting (64)" errBits2;
                             {0}
                             '""".format(extracmd),
                 'out': output_filename}

        cmd = '{stilts} tcat {in} icmd={icmd} countrows=true lazy=true out={out}'
        mycmd = cmd.format(**param)
        log.info(mycmd)
        status = os.system(mycmd)
        log.info('concat: '+str(status))

        # zip
        mycmd = 'gzip --stdout {0} > {1}'.format(output_filename, output_filename_gzip)
        log.debug(mycmd)
        status = os.system(mycmd)
        log.info('gzip: '+str(status))

        return status


###########
# FUNCTIONS
###########

def concatenate_one(strip,
                    logfile = os.path.join(constants.LOGDIR, 'concatenation.log')):
    with log.log_to_file(logfile):
        # Strips are defined by the start longitude
        log.info('Concatenating L={0}'.format(strip))
        for mode in ['light', 'full']:
            for part in ['a', 'b']:
                concat = Concatenator(strip, part, mode)
                concat.run()
    return strip


def concatenate(clusterview):
    # Spread the work across the cluster
    strips = np.arange(25, 215+1, constants.STRIPWIDTH)
    results = clusterview.imap(concatenate_one, strips)

    # Print a friendly message once in a while
    i = 0
    for mystrip in results:
        i += 1
        log.info('Completed strip {0} ({1}/{2})'.format(mystrip,
                                                        i,
                                                        len(strips)))
    log.info('Concatenating finished')


def merge_light_catalogue():
    """Merge the light tiled catalogues into one big file."""
    output_filename = os.path.join(constants.DESTINATION,
                                   'concatenated',
                                   'iphas-dr2-light.fits')

    instring = ''
    for lon in np.arange(25, 215+1, constants.STRIPWIDTH):
        for part in ['a', 'b']:
            path = os.path.join(constants.DESTINATION,
                                'concatenated',
                                'light',
                                'iphas-dr2-{0:03d}{1}-light.fits'.format(
                                                                    lon, part))
            instring += 'in={0} '.format(path)

    # Warning: a bug in stilts causes long fieldIDs to be truncated if -utype S15 is not set
    param = {'stilts': constants.STILTS,
             'in': instring,
             'out': output_filename}

    cmd = '{stilts} tcat {in} countrows=true lazy=true ofmt=colfits-basic out={out}'
    mycmd = cmd.format(**param)
    log.debug(mycmd)
    status = os.system(mycmd)
    log.info('concat: '+str(status))

    return status


################################
# MAIN EXECUTION (FOR DEBUGGING)
################################

if __name__ == "__main__":
    log.setLevel('DEBUG')
    concatenate_one(215)
