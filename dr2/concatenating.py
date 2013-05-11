#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Concatenates the source lists into the final Source Catalogue.

This script takes the output from the seaming script and concatenates the
results into the final "Primary Source Catalogue" product,
which is generated in 5x5 degree tiles.

Both a light and a full version are generated.

"""
from __future__ import division, print_function, unicode_literals
import os
import numpy as np
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

    def __init__(self, strip, part='a', mode='full', calibrated=True):
        assert(part in ['a', 'b'])
        assert(mode in ['light', 'full'])
        assert(calibrated in [True, False])

        self.strip = strip
        self.part = part
        self.mode = mode
        self.calibrated = calibrated

        # Where are the input catalogues?
        if calibrated:
            self.datapath = os.path.join(constants.DESTINATION, 'calibrated')
        else:
            self.datapath = os.path.join(constants.DESTINATION, 'seamed')

        # Where to write the output?
        if self.calibrated:
            self.destination = os.path.join(constants.DESTINATION,
                                            'concatenated')
        else:
            self.destination = os.path.join(constants.DESTINATION,
                                            'concatenated-uncalibrated')
        if mode == 'light':
            self.destination = os.path.join(self.destination, 'light')

        # Make sure our destination exists
        if not os.path.exists(self.destination):
            os.makedirs(self.destination)

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
            suffix = '-full'

        return os.path.join(self.destination,
                            'iphas-dr2-psc-glon{0:03.0f}{1}{2}.fits'.format(
                                                    self.lon1,
                                                    self.part,
                                                    suffix))

    def get_fieldlist(self):
        # Which are our fields?
        # Note: we must allow for border overlaps
        if self.part == 'a':
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
        if self.part == 'a':
            cond_latitude = "b >= 0"
        else:
            cond_latitude = "b < 0"

        if self.mode == 'full':
            extracmd = """delcols "rPlaneX rPlaneY iPlaneX iPlaneY \
                                   haPlaneX haPlaneY rAxis primaryID" """
        else:
            extracmd = """select "nBands == 3"; \
                          keepcols "ra dec \
                                       r rErr \
                                       i iErr \
                                       ha haErr \
                                       mergedClass errBits";"""

        instring = ''
        for field in self.fieldlist:
            path = os.path.join(self.datapath,
                                'strip{0}'.format(self.strip),
                                '{0}.fits'.format(field))
            instring += 'in={0} '.format(path)

        output_filename = self.get_output_filename()
        log.info('Writing data to {0}'.format(output_filename))

        # A bug in stilts causes long fieldIDs to be truncated if -utype S15 is not set
        param = {'stilts': constants.STILTS,
                 'in': instring,
                 'icmd': """'select "(errBits < 64) \
                                      & (rErr < 0.198 || iErr < 0.198 || haErr < 0.198)
                                      & (pStar > 0.2 || pGalaxy > 0.2) \
                                      & l >= """+str(self.lon1)+""" \
                                      & l < """+str(self.lon2)+""" \
                                      & """+str(cond_latitude)+""" \
                                      & sourceID == primaryID"; \
                             replacecol -utype S15 fieldID "fieldID"; \
                             replacecol -utype S1 fieldGrade "toString(fieldGrade)"; \
                             {0}
                             '""".format(extracmd),
                 'out': output_filename}

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
    for mode in ['light', 'full']:
        for part in ['a', 'b']:
            concat = Concatenator(strip, part, mode, calibrated=True)
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
        run_all(2)
