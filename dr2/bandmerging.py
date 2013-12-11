#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generates the user-friendly band-merged catalogues.

This script will merge the same-epoch H-alpha/r/i exposures of each IPHAS
field into band-merged catalogues. The column names follow the UKIDSS conventions.
"""
from __future__ import print_function, division
import os
import sys
import numpy as np
from astropy import log
from multiprocessing import Pool
import constants
import util
from constants import IPHASQC
import socket

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Geert Barentsen', 'Hywel Farnhill', 
               'Janet Drew', 'Robert Greimel']


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
        self.output = os.path.join(MYDESTINATION, fieldid.strip()+'.fits')

    def get_catalogue_path(self, run):
        """Returns the full path of a run's 'iphasDetection' catalogue."""
        return os.path.join(constants.DESTINATION,
                            'detected',
                            '{}_det.fits'.format(run))

    def get_stilts_command(self):
        """Returns the stilts command used to perform a band-merge."""

        config = {'stilts': constants.STILTS,
                  'MATCHING_DISTANCE': constants.MATCHING_DISTANCE,
                  'runr': self.get_catalogue_path(self.run_r),
                  'runi': self.get_catalogue_path(self.run_i),
                  'runha': self.get_catalogue_path(self.run_ha),
                  'fieldgrade': self.fieldgrade,
                  'fieldid': self.fieldid,
                  'ocmd': os.path.join(constants.LIBDIR,
                                       'stilts-band-merging.cmd'),
                  'output': self.output}

        # Note: we filter out sources detected at <0.5sigma
        # (i.e. magnitude error < -2.5*log(1+3) = 1.19)
        cmd = """{stilts} tmatchn matcher=sky params={MATCHING_DISTANCE} \
                  multimode=group nin=3 \
                  in1={runr} in2={runi} in3={runha} \
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
        if status != 0:
            log.error('{0}: Unexpected status ("{1}"): command was: {2}'.format(self.fieldid, status, cmd))
        return status


###########
# FUNCTIONS
###########

def bandmerge_one(fieldid):
    """Band-merge a single field """
    with log.log_to_file(os.path.join(constants.LOGDIR, 'bandmerging.log')):
        engine = socket.gethostname()+'/'+str(os.getpid())
        log.info('Starting {0} on {1}'.format(fieldid, engine))

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

        log.info('Finished {0} on {1} (returned {2})'.format(fieldid,
                                                             engine,
                                                             status))
        return status


def bandmerge(clusterview):
    """Band-merge all fields."""
    util.setup_dir(MYDESTINATION)

    # Spread the work across the cluster
    field_ids = IPHASQC.field('id')
    results = clusterview.imap(bandmerge_one, field_ids)

    # Print a friendly message once in a while
    i = 0
    for status in results:
        i += 1
        if (i % 1000) == 0:
            log.info('Completed field {0}/{1}'.format(i, len(field_ids)))
    log.info('Bandmerging finished')



###################
# MAIN EXECUTION
###################

if __name__ == '__main__':
    if constants.DEBUGMODE:
        log.setLevel('INFO')
        #run_all(lon1=lon1, lon2=lon2, ncores=4)
        #bandmerge_one('5089o_jun2005')
        #bandmerge_one('3561_nov2003')
