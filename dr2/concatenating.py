#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Concatenates the source lists into the final Source Catalogue.

This script takes the output from the seaming script and concatenates the
results into the final "Primary Source Catalogue" product,
which is generated in 5x5 degree tiles.

TODO
 * Offer a 'full' and 'light' version?
"""
from __future__ import division, print_function, unicode_literals
import os
import numpy as np
from multiprocessing import Pool
from astropy.io import fits
from astropy import log
import constants
from constants import IPHASQC

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Hywel Farnhill', 'Robert Greimel', 'Janet Drew']


#############################
# CONSTANTS & CONFIGURATION
#############################

# Where to write the output catalogues?
MYDESTINATION = os.path.join(constants.DESTINATION, 'concatenated')


###########
# CLASSES
###########

class Concatenator(object):

    def __init__(self, strip, mode):
        self.strip = strip
        self.mode = mode
        assert(self.mode in ['a', 'b'])
        self.lon1 = strip
        self.lon2 = strip + constants.STRIPWIDTH
        self.fieldlist = self.get_fieldlist()

        if not os.path.exists(MYDESTINATION):
            os.makedirs(MYDESTINATION)

        self.output_file = os.path.join(MYDESTINATION,
                                        'iphas-dr2-psc-glon{0:03.0f}{1}.fits'.format(
                                                    self.lon1, mode))

    def get_fieldlist(self):
        # Which are our fields?
        # Note: we must allow for border overlaps
        if self.mode == 'a':
            cond_b = IPHASQC['b'] > (0 - constants.FIELD_MAXDIST)
        else:
            cond_b = IPHASQC['b'] < (0 + constants.FIELD_MAXDIST)

        cond_strip = (constants.IPHASQC_COND_RELEASE
                      & cond_b
                      & (IPHASQC['l'] >= (self.lon1 - constants.FIELD_MAXDIST))
                      & (IPHASQC['l'] < (self.lon2 + constants.FIELD_MAXDIST)))
        fieldlist = IPHASQC['id'][cond_strip]
        log.info('Found {0} fields.'.format(len(fieldlist)))
        return fieldlist

    def run(self):
        """Performs the concatenation of the strip."""
        if self.mode == 'a':
            cond_latitude = "b >= 0"
        else:
            cond_latitude = "b < 0"

        instring = ''
        for field in self.fieldlist:
            path = os.path.join(constants.DESTINATION,
                                'seamed',
                                'strip{0}'.format(self.strip),
                                '{0}.fits'.format(field))
            instring += 'in={0} '.format(path)
        # & (sourceID == primaryID) \
        # Note: a bug in stilts causes long fieldIDs to be truncated if -utype S15 is not set
        param = {'stilts': constants.STILTS,
                 'in': instring,
                 'icmd': """'select "(errBits < 100) \
                                      & (pStar > 0.2 || pGalaxy > 0.2) \
                                      & l >= """+str(self.lon1)+""" \
                                      & l < """+str(self.lon2)+""" \
                                      & """+str(cond_latitude)+""" \
                                      & sourceID == primaryID"; \
                             replacecol -utype S15 fieldID "fieldID"; \
                             replacecol -utype S1 fieldGrade "toString(fieldGrade)"; \
                             keepcols "sourceID ra dec l b \
                                       mergedClass mergedClassStat \
                                       pStar pGalaxy \
                                       rmi rmha \
                                       r rErr rAperMag3 rClass rMJD rErrBits \
                                       i iErr iAperMag3 iClass iMJD iErrBits iXi iEta \
                                       ha haErr haAperMag3 haClass haMJD \
                                       haXi haEta \
                                       haPeakMag haPeakMagErr haClassStat \
                                       brightNeighb deblend saturated \
                                       vignetted errBits reliable reliableStar \
                                       fieldID fieldGrade night seeing \
                                       nObs r2 rErr2 i2 iErr2 \
                                       ha2 haErr2 errBits2"'""",
                 'out': self.output_file}

        cmd = '{stilts} tcat {in} icmd={icmd} countrows=true lazy=true out={out}'
        mycmd = cmd.format(**param)
        log.debug(mycmd)
        status = os.system(mycmd)
        log.info('Concat: '+str(status))
        return status




###########
# FUNCTIONS
###########

def run_one(strip):
    # Strips are defined by the start longitude of a 10 deg-wide strip
    #assert(strip in np.arange(30, 210+1, 10))
    log.info('Concatenating L={0}'.format(strip))
    concat = Concatenator(strip, 'a')  # north
    concat.run()
    concat = Concatenator(strip, 'b')  # south
    concat.run()


def run_all(ncores=2):
    longitudes = np.arange(25, 215+1, constants.STRIPWIDTH)[::-1]
    # Run the processing for each pipeline catalogue
    p = Pool(processes=ncores)
    p.map(run_one, longitudes)

###################
# MAIN EXECUTION
###################

if __name__ == "__main__":
    if constants.DEBUGMODE:
        log.setLevel('INFO')
        run_one(200)
    else:
        log.setLevel('INFO')
        run_all(4)
