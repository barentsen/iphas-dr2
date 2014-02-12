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
            extracmd = """delcols "rPlaneX rPlaneY iPlaneX iPlaneY \
                                   haPlaneX haPlaneY rAxis primaryID
                                   vignetted truncated badPix" """
        else:
            # select "nBands == 3"; \
            extracmd = """replacecol errBits "toShort(errBits)";
                          keepcols "name ra dec \
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
                                    -desc "Source designation, JHHMMSS.ss+DDMMSS.s" \
                                    name \
                                    "concat(\\"J\\", 
                                            replaceAll(degreesToHms(ra, 2),
                                                       \\":\\", \\"\\"), 
                                            replaceAll(degreesToDms(dec, 1),
                                                       \\":\\", \\"\\")
                                            )"; \
                             replacecol -utype S15 fieldID "fieldID"; \
                             replacecol -utype S1 fieldGrade "toString(fieldGrade)"; \
                             replacecol errBits "toShort(errBits)";
                             replacecol rErrBits "toShort(rErrBits)";
                             replacecol iErrBits "toShort(iErrBits)";
                             replacecol haErrBits "toShort(haErrBits)";
                             colmeta -desc "Human-readable IPHAS field number and observing run (e.g. 0001o_aug2003)." fieldID;
                             colmeta -desc "Internal quality control score of the field. One of A, B, C or D." fieldGrade;
                             colmeta -desc "Number of repeat observations of this source in the survey." nObs;
                             colmeta -desc "SourceID of the object in the partner exposure (if obtained within 10 minutes of the primary detection)." sourceID2;
                             colmeta -desc "FieldID of the partner detection (e.g. 0001o_aug2003)." fieldID2;
                             colmeta -desc "r-band magnitude in the partner field, i.e. the dithered repeat measurement obtained within 10 minutes (if available)." r2;
                             colmeta -desc "Uncertainty for r2." rErr2;
                             colmeta -desc "i-band magnitude in the partner field, i.e. the dithered repeat measurement obtained within 10 minutes (if available)." i2;
                             colmeta -desc "Uncertainty for i2." iErr2;
                             colmeta -desc "H-alpha magnitude in the dithered partner field, i.e. the dithered repeat measurement obtained within 10 minutes (if available)." ha2;
                             colmeta -desc "Uncertainty for ha2." haErr2;
                             colmeta -desc "Error bitmask for the partner detection. Used to flag a bright neighbour (1), source blending (2), saturation (8), vignetting (64), truncation (128) and bad pixels (32768)." errBits2;
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


###################
# MAIN EXECUTION
###################

if __name__ == "__main__":
    log.setLevel('DEBUG')
    concatenate_one(215)
