#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generates the user-friendly band-merged 'iphasSource' catalogues.

This script will merge the same-epoch H-alpha/r/i exposures of each IPHAS
field into band-merged catalogues, following the UKIDSS column definitions.

TODO:
 * add htmID?
"""
from __future__ import print_function, division
import os
import sys
import numpy as np
from astropy import log
from multiprocessing import Pool
import constants
from constants import IPHASQC

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Hywel Farnhill', 'Robert Greimel', 'Janet Drew']


#############################
# CONSTANTS & CONFIGURATION
#############################

# Where to write the output catalogues?
MYDESTINATION = os.path.join(constants.DESTINATION, 'bandmerged')


###########
# CLASSES
###########

class IPHASException(Exception):
    """Raised when a field cannot be band-merged."""
    pass


class BandMerge():
    """
    Class to read in iphasSource tables and perform a band-merge.
    """

    def __init__(self, fieldid, fieldgrade, run_ha, run_r, run_i):
        """Constructor"""
        self.fieldid = fieldid
        self.fieldgrade = fieldgrade
        self.run_ha = run_ha
        self.run_r = run_r
        self.run_i = run_i
        self.output = os.path.join(MYDESTINATION, fieldid+'.fits')

    def get_catalogue_path(self, run):
        """Returns the full path of a run's 'iphasDetection' catalogue."""
        return os.path.join(constants.DESTINATION,
                            'detected',
                            '{}_det.fits'.format(run))

    def get_stilts_command(self):
        """Returns the stilts command used to perform a band-merge."""

        config = {'stilts': constants.STILTS,
                  'runr': self.get_catalogue_path(self.run_r),
                  'runi': self.get_catalogue_path(self.run_i),
                  'runha': self.get_catalogue_path(self.run_ha),
                  'fieldgrade': self.fieldgrade,
                  'fieldid': self.fieldid,
                  'ocmd': os.path.join(constants.PACKAGEDIR,
                                       'lib',
                                       'stilts-band-merging.cmd'),
                  'output': self.output}

        # Note: we filter out sources detected at <0.5sigma
        # (i.e. magnitude error < -2.5*log(1+3) = 1.19)
        cmd = """{stilts} tmatchn matcher=sky params=0.5 multimode=group \
                  nin=3 in1={runr} in2={runi} in3={runha} \
                  join1=always join2=always join3=always \
                  values1='ra dec' values2='ra dec' values3='ra dec' \
                  icmd1='setparam fieldID "{fieldid}";
                         setparam fieldGrade "{fieldgrade}";
                         select "aperMag2Err > 0 & aperMag2Err < 1.19";' \
                  icmd2='select "aperMag2Err > 0 & aperMag2Err < 1.19";' \
                  icmd3='select "aperMag2Err > 0 & aperMag2Err < 1.19";' \
                  ocmd=@{ocmd} \
                  progress=none \
                  out='{output}'""".format(**config)

        return cmd

    def run(self):
        """Perform the band-merging."""
        cmd = self.get_stilts_command()
        status = os.system(cmd)
        return status


###########
# FUNCTIONS
###########

def run_one(fieldid):
    """ Band-merge a single field """
    # Which index does the field have in the QC table?
    idx = np.where(IPHASQC.field('id') == fieldid)[0]
    if len(idx) < 1:
        raise IPHASException('{}: error identifying runs'.format(fieldid))

    # Carry out the band-merging
    bm = BandMerge(fieldid,
                   IPHASQC.field('qflag')[idx[0]],
                   IPHASQC.field('run_ha')[idx[0]],
                   IPHASQC.field('run_r')[idx[0]],
                   IPHASQC.field('run_i')[idx[0]])
    status = bm.run()

    log.info('{0}: {1}'.format(fieldid, status))
    return status


def run_all(lon1=20, lon2=220, ncores=4):
    """ Band-merge all fields """
    log.info('Band-merging in longitude range [{0},{1}['.format(lon1, lon2))

    # Make sure the output directory exists
    if not os.path.exists(MYDESTINATION):
        os.makedirs(MYDESTINATION)

    # Which fields do we want to merge?
    idx = np.where(constants.IPHASQC_COND_RELEASE
                   & (IPHASQC.field('l') >= lon1)
                   & (IPHASQC.field('l') < lon2))[0]
    fieldids = IPHASQC.field('id')[idx]
    # Distribute the work over ncores
    p = Pool(processes=ncores)
    results = p.imap(run_one, fieldids)
    for i in results:
        pass


###################
# MAIN EXECUTION
###################

if __name__ == '__main__':
    """Bandmerges fields in the IPHAS survey within a given gal. long. range.

    Arguments
    ---------
    lon1: number, default 0
    begin galactic longitude (default: 0)
    lon2: number, default 360
    end galactic longitude (default: 360)
    """

    # Which longitude strip to process?
    if len(sys.argv) > 2:
        lon1 = int(sys.argv[1])
        lon2 = int(sys.argv[2])
    else:  # Do all
        lon1 = 0
        lon2 = 360

    if constants.DEBUGMODE:
        log.setLevel('INFO')
        #run_all(lon1=lon1, lon2=lon2, ncores=4)
        #run_one('5089o_jun2005')
        #run_one('3561_nov2003')
        #run_all(lon1=208, lon2=209, ncores=6)

        run_one('2925o_dec2003')
        run_one('2832_dec2007')
        run_one('2914o_dec2003')
        run_one('2327_dec2007')
        run_one('2426_oct2006')
        run_one('6745_sep2005')
        run_one('2798_nov2012')
        run_one('3367o_oct2006')
        run_one('3352o_oct2006')
        run_one('1385o_oct2004')
        run_one('2694o_dec2005')
        run_one('2817o_dec2003')
        run_one('2694o_dec2005')
        run_one('2831o_dec2003')

    else:  # production
        run_all(lon1=lon1, lon2=lon2, ncores=8)
