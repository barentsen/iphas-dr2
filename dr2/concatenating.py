#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Concatenates the bandmerged field catalogues into a "Point Source Catalogue".

This script takes the output from the seaming script and concatenates the
results into the final "Primary Source Catalogue" product,
which is generated in 5x5 degree tiles.

Both a light and a full version are produced.
"""
from __future__ import division, print_function, unicode_literals
import os
import numpy as np
import datetime
from multiprocessing import Pool
from astropy import log
import constants
from constants import IPHASQC

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Hywel Farnhill', 'Robert Greimel', 'Janet Drew']


###########
# CLASSES
###########

class Concatenator(object):

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
        
        if mode == 'light':
            self.destination = os.path.join(self.destination, 'light')
        else:
            self.destination = os.path.join(self.destination, 'full')

        # Make sure our destination exists
        try:
            if not os.path.exists(self.destination):
                os.makedirs(self.destination)
        except Exception:
            pass

        log.info('Reading data from {0}'.format(self.datapath))

        # Limits
        self.lon1 = strip
        self.lon2 = strip + constants.STRIPWIDTH
        self.fieldlist = self.get_fieldlist()

    def get_output_filename(self):
        """Returns the full path of the output file."""
        if self.mode == 'light':
            suffix = '-light'
        else:
            suffix = ''

        return os.path.join(self.destination,
                            'iphas-dr2-{0:03.0f}{1}{2}.fits'.format(
                                                    self.lon1,
                                                    self.part,
                                                    suffix))

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
        """Performs the concatenation of the strip."""
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
            extracmd = """replacecol ra "toFloat(ra)";
                          replacecol dec "toFloat(dec)";
                          replacecol errBits "toShort(errBits)";
                          keepcols "ra dec \
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
        log.info('Writing data to {0}'.format(output_filename))

        version = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')

        # A bug in stilts causes long fieldIDs to be truncated if -utype S15 is not set
        param = {'stilts': constants.STILTS,
                 'in': instring,
                 'icmd': """'clearparams *; \
                             setparam NAME "IPHAS DR2 Source Catalogue"; \
                             setparam ORIGIN "www.iphas.org"; \
                             setparam AUTHOR "Geert Barentsen, Hywel Farnhill, Janet Drew"; \
                             setparam VERSION \""""+version+""""; \
                             select "(errBits < 64) \
                                      & (rErr < 0.198 || iErr < 0.198 || haErr < 0.198) \
                                      & (pStar > 0.2 || pGalaxy > 0.2) \
                                      & l >= """+str(self.lon1)+""" \
                                      & l < """+str(self.lon2)+""" \
                                      & """+str(cond_latitude)+""" \
                                      & sourceID == primaryID"; \
                             addcol -before ra \
                                    -desc "Source designation, JHHMMSS.ss+DDMMSS.s" \
                                    name \
                                    "concat(\"J\", replaceAll(degreesToHms(ra, 2), \":\", \"\"), replaceAll(degreesToDms(dec, 1), \":\", \"\"))"; \
                             replacecol -utype S15 fieldID "fieldID"; \
                             replacecol -utype S1 fieldGrade "toString(fieldGrade)"; \
                             {0}
                             '""".format(extracmd),
                 'out': output_filename}

        cmd = '{stilts} tcat {in} icmd={icmd} countrows=true lazy=true out={out}'
        mycmd = cmd.format(**param)
        log.info(mycmd)
        status = os.system(mycmd)
        log.info('concat: '+str(status))

        # zip
        mycmd = 'gzip {0}'.format(output_filename)
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
