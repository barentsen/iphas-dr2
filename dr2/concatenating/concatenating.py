#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generates the final IPHAS data products."""
from __future__ import division, print_function, unicode_literals
import os
import sys
import numpy as np
from astropy.io import fits
from astropy import log
log.setLevel('DEBUG')

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

    def __init__(self, fieldlist, strip):
        self.fieldlist = fieldlist
        self.strip = strip

        self.lon1 = strip
        self.lon2 = strip + 10
        if self.lon1 == 30:
            self.lon1 = 25

        if not os.path.exists(DESTINATION):
            os.makedirs(DESTINATION)

        self.output_file = os.path.join(DESTINATION,
                                        'strip{0}.fits'.format(strip))

    def run(self):
        """Performs the concatenation of the strip."""
        instring = ''
        for field in self.fieldlist:
            path = os.path.join(DATADIR, 
                               'strip{0}'.format(self.strip),
                               '{0}.fits'.format(field))
            instring += 'in={0} '.format(path)
        # & (sourceID == primaryID) \
        param = {'stilts': STILTS,
                 'in': instring,
                 'icmd': """'select "(errBits < 100) \
                                      & pStar > 0.2 \
                                      & l >= """+str(self.lon1)+""" \
                                      & l < """+str(self.lon2)+""""; \
                             keepcols "sourceID primaryID ra dec l b mergedClass \
                                       r rErr rClass i iErr iClass ha haErr haClass \
                                       brightNeighb deblend saturated vignetted \
                                       errBits reliable reliableStar fieldID"'""",
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

def run_strip(strip):
    # Strips are defined by the start longitude of a 10 deg-wide strip
    assert(strip in np.arange(30, 210+1, 10))
    log.info('Concatenating strip {0}'.format(strip))

    # So which are our boundaries?
    # Note: we must allow an extree for border overlaps!
    lon1 = strip - 0.8
    lon2 = strip + 10 + 0.8
    cond_strip = (IPHASQC['is_pdr']
                  & (IPHASQC['l'] >= lon1)
                  & (IPHASQC['l'] < lon2))

    fieldlist = IPHASQC['id'][cond_strip]
    log.info('Found {0} fields.'.format(len(fieldlist)))

    concat = Concatenator(fieldlist, strip)
    concat.run()



###################
# MAIN EXECUTION
###################

if __name__ == "__main__":

    # Which longitude range to process?
    if len(sys.argv) > 1:
        strip = int(sys.argv[1])
    else:
        raise Exception('Missing longitude strip argument')

    run_strip(strip)
