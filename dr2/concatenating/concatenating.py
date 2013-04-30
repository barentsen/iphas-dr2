#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generates the final IPHAS data products.

TODO:
 * fieldID appears truncated if the first row has a short fieldID;
"""
from __future__ import division, print_function, unicode_literals
import os
import sys
import numpy as np
from multiprocessing import Pool
from astropy.io import fits
from astropy import log
log.setLevel('INFO')

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Hywel Farnhill', 'Robert Greimel', 'Janet Drew',
               'Cambridge Astronomical Survey Unit']


#############################
# CONSTANTS & CONFIGURATION
#############################

HOSTNAME = os.uname()[1]
if HOSTNAME == 'uhppc11.herts.ac.uk':  # testing
    # Where are the pipeline-reduced catalogues?
    DATADIR = "/home/gb/tmp/iphas-dr2/seamed"
    # Where to write the output catalogues?
    DESTINATION = "/home/gb/tmp/iphas-dr2/concatenated"
else:  # production
    DATADIR = "/car-data/gb/iphas-dr2/seamed"
    DESTINATION = "/car-data/gb/iphas-dr2/concatenated"

# Where is the IPHAS quality control table?
IPHASQC = fits.getdata('/home/gb/dev/iphas-qc/qcdata/iphas-qc.fits', 1)

# Dir of this script
SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))

# How to execute stilts?
STILTS = 'nice java -Xmx2000M -XX:+UseConcMarkSweepGC -jar {0}'.format(
                                 os.path.join(SCRIPTDIR, '../lib/stilts.jar'))

###########
# CLASSES
###########

class Concatenator(object):

    def __init__(self, strip, lon1, lon2):
        self.strip = strip
        self.lon1 = lon1
        self.lon2 = lon2
        self.fieldlist = self.get_fieldlist()

        if not os.path.exists(DESTINATION):
            os.makedirs(DESTINATION)

        self.output_file = os.path.join(DESTINATION,
                                        'iphas-dr2-psc-glon{0:03.0f}.fits'.format(
                                                    self.lon1))

    def get_fieldlist(self):
        # Which are our fields?
        # Note: we must allow for border overlaps
        cond_strip = (IPHASQC['is_pdr']
                      & (IPHASQC['qflag'] != 'D')
                      & (IPHASQC['l'] >= (self.lon1 - 0.8))
                      & (IPHASQC['l'] < (self.lon2 + 0.8)))
        fieldlist = IPHASQC['id'][cond_strip]
        log.info('Found {0} fields.'.format(len(fieldlist)))
        return fieldlist

    def run(self):
        """Performs the concatenation of the strip."""
        instring = ''
        for field in self.fieldlist:
            path = os.path.join(DATADIR,
                               'strip{0}'.format(self.strip),
                               '{0}.fits'.format(field))
            instring += 'in={0} '.format(path)
        # & (sourceID == primaryID) \
        # Note: a bug in stilts causes long fieldIDs to be truncated if -utype S15 is not set
        param = {'stilts': STILTS,
                 'in': instring,
                 'icmd': """'select "(errBits < 100) \
                                      & (pStar > 0.2 || pGalaxy > 0.2) \
                                      & l >= """+str(self.lon1)+""" \
                                      & l < """+str(self.lon2)+"""
                                      & sourceID == primaryID"; \
                             replacecol -utype S15 fieldID "fieldID"; \
                             replacecol -utype S1 fieldGrade "toString(fieldGrade)"; \
                             keepcols "sourceID ra dec l b \
                                       mergedClass mergedClassStat \
                                       pStar pGalaxy \
                                       rmi rmha \
                                       r rErr rAperMag3 rClass rMJD \
                                       i iErr iAperMag3 iClass iMJD iXi iEta \
                                       ha haErr haAperMag3 haClass haMJD haXi haEta \
                                       haPeakMag haPeakMagErr haClassStat \
                                       brightNeighb deblend saturated vignetted \
                                       errBits reliable reliableStar \
                                       fieldID fieldGrade night seeing"'""",
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

def run_one(lon1):
    if lon1 < 30:
        strip = 30
    else:
        strip = lon1 - (lon1 % 10)

    lon2 = lon1 + 5

    # Strips are defined by the start longitude of a 10 deg-wide strip
    #assert(strip in np.arange(30, 210+1, 10))
    log.info('Concatenating L={0}-{1}'.format(lon1, lon2))
    concat = Concatenator(strip, lon1, lon2)
    concat.run()


def run_all(ncores=8):
    longitudes = np.arange(25, 215+1, 5)[::-1]
    # Run the processing for each pipeline catalogue
    p = Pool(processes=ncores)
    p.map(run_one, longitudes)

###################
# MAIN EXECUTION
###################

if __name__ == "__main__":

    run_all(8)
    #run_one(215)
