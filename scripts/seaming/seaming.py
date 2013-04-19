#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Find the primary detections in the IPHAS catalogue."""
from __future__ import division, print_function, unicode_literals
import os
import sys
import numpy as np
from astropy.io import fits
import logging
logging.basicConfig(level=logging.INFO)

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Hywel Farnhill', 'Robert Greimel', 'Janet Drew',
               'Cambridge Astronomical Survey Unit']


#############################
# CONSTANTS & CONFIGURATION
#############################

# Where is the IPHAS quality control table?
IPHASQC = fits.getdata('/home/gb/dev/iphas-qc/qcdata/iphas-qc.fits', 1)

# How to execute stilts?
STILTS = 'nice java -Xmx500M -XX:+UseConcMarkSweepGC -jar ../lib/stilts.jar'

FIELD_MAXDIST = 0.8  # Fields within this radius will be considered to overlap

CACHE = {}


###########
# CLASSES
###########

class SeamingException(Exception):
    """
    Exception raised when a catalogue has a known problem.
    """
    pass


class Seamer(object):

    def __init__(self, fieldid, ra, dec):
        self.fieldid = fieldid
        self.crossmatch_file = '/tmp/seaming_{0}.fits'.format(fieldid)
        self.primaryid_file = '/tmp/primaryid_{0}.fits'.format(fieldid)
        self.output_file = '/home/gb/tmp/iphas-dr2/iphasSource2/{0}.fits'.format(fieldid)
        self.ra = ra
        self.dec = dec

    def run(self):
        self.overlaps = self.get_overlaps()
        self.crossmatch()
        sourceID, primaryID = self.get_primary_ids()
        #print(primaryID)
        self.add_column(sourceID, primaryID)
        self.clean()


    def clean(self):
        os.remove(self.crossmatch_file)
        os.remove(self.primaryid_file)


    def add_column(self, sourceID, primaryID):
        # Write the (sourceID,primaryID)s to a table
        col_sourceID = fits.Column(name='sourceID', format='K', array=sourceID)
        col_priSourceID = fits.Column(name='priSourceID', format='K', array=primaryID)
        cols = fits.ColDefs([col_sourceID, col_priSourceID])
        newtable = fits.new_table(cols)
        newtable.writeto(self.primaryid_file, clobber=True)

        # Then use stilts to add the extra column
        config = {'STILTS': STILTS,
                  'IN1': self.filename(self.fieldid),
                  'IN2': self.primaryid_file,
                  'OUT': self.output_file}

        cmd = "{STILTS} tmatch2 progress=none find=best1 in1={IN1} in2={IN2} "
        cmd += "matcher=exact join=all1 suffix1='' "
        cmd += "values1='sourceID' values2='sourceID' "
        cmd += "ocmd='delcols sourceID_2' out='{OUT}' "

        stilts_cmd = cmd.format(**config)
        logging.info(stilts_cmd)
        status = os.system(stilts_cmd)
        logging.info(status)
        return status

    def get_overlaps(self):
        """Which fields possibly overlap?"""
        ddec = self.dec - IPHASQC['dec']
        dra = (self.ra - IPHASQC['ra']) * np.cos(np.radians(self.dec))
        d = (ddec ** 2 + dra ** 2) ** 0.5
        idx = (IPHASQC['is_pdr']
               & (d < FIELD_MAXDIST)
               & (self.fieldid != IPHASQC['id']))
        return IPHASQC['id'][idx]

    def filename(self, fieldid):
        """
        url = ('http://stri-cluster.herts.ac.uk/~gb/'
               + 'iphas-dr2/iphasSource/{0}.fits')
        """
        url = '/home/gb/tmp/iphas-dr2/iphasSource/{0}.fits'
        return url.format(fieldid)

    def get_stilts_command(self):
        # Operations to perform on all tables
        addcol = 'addcol bands "sum(array(NULL_r?0:1,NULL_i?0:1,NULL_ha?0:1))"'
        keepcols = 'keepcols "sourceID fieldID ra dec bands rSeeing"'
        # Keywords in stilts command
        config = {'STILTS': STILTS,
                  'NIN': len(self.overlaps) + 1,
                  'IN1': self.filename(self.fieldid),
                  'ICMD': addcol+'; '+keepcols,
                  'OUT': self.crossmatch_file}
        # Create stilts command
        # FIXME: add progress=none
        cmd = "{STILTS} tmatchn progress=none matcher=sky params=0.5 "
        cmd += "multimode=pairs iref=1 nin={NIN} "
        cmd += "in1={IN1} values1='ra dec' join1=always "
        cmd += "icmd1='{ICMD}' "
        for i in range(len(self.overlaps)):
            cmd += "in{0}={1} ".format(i+2,
                                       self.filename(self.overlaps[i]))
            cmd += "values{0}='ra dec' ".format(i+2)
            cmd += "icmd{0}='{ICMD}' ".format(i+2, **config)
        cmd += "out='{OUT}' "

        stilts_cmd = cmd.format(**config)
        return stilts_cmd

    def crossmatch(self):
        cmd = self.get_stilts_command()
        logging.debug('Overlaps: {0}'.format(self.overlaps))
        logging.debug(cmd)
        status = os.system(cmd)
        logging.info(status)
        return status

    def get_primary_ids(self):
        try:
            d = fits.getdata(self.crossmatch_file, 1)
        except IOError, e:
            raise SeamingException(e)

        primary_ids = []
        for source in d:

            if source['sourceID_1'] in CACHE:
                 primary_ids.append(CACHE[source['sourceID_1']])
            else:

                idx_all = (np.arange(self.overlaps.size + 1) + 1)
                ids = np.array([source['sourceID_{0}'.format(i)] for i in idx_all])

                idx = idx_all[ids > 0]

                stats = {'sourceID': ids[ids > 0],
                         'fieldID': np.array([source['fieldID_{0}'.format(i)]
                                              for i in idx]),
                         'bands': np.array([source['bands_{0}'.format(i)]
                                            for i in idx]),
                         'errBits': np.array([0 for i in idx]),  # FIXME
                         'seeing': np.array([source['rSeeing_{0}'.format(i)]
                                             for i in idx]),  # FIXME
                         'rAxis': np.array([i for i in idx]),  # FIXME
                         'winner': np.array([True for i in idx])}

                stats = self.mark_winner(stats)

                winnerID = stats['sourceID'][stats['winner']][0]

                primary_ids.append(winnerID)

                for sourceID in stats['sourceID']:
                    CACHE[sourceID] = winnerID

        return (d['sourceID_1'], np.array(primary_ids))

    def mark_winner(self, stats):
        """The ranking criteria are
         1. filter coverage (3 > 2 > 1);
         2. error flags;
         3. worst seeing amongst the filters 
            (keep all detections within 20% of the best value);
         4. distance from the detector+ccd edge.
        """
        if stats['winner'].sum() > 1:

            # Discard source with less bands
            idx_loser = (stats['bands'] < np.max(stats['bands']))
            stats['winner'][idx_loser] = False

            # Discard sources with large errBits
            idx_loser = (stats['errBits']
                         > np.min(stats['errBits'][stats['winner']]))
            stats['winner'][idx_loser] = False

            # Discard sources with poor seeing
            idx_loser = (stats['seeing']
                         > 1.2 * np.min(stats['seeing'][stats['winner']]))
            stats['winner'][idx_loser] = False

            # Finally, the winner is the one closest to the optical axis
            idx_loser = (stats['rAxis']
                         > np.min(stats['rAxis'][stats['winner']]))
            stats['winner'][idx_loser] = False

            if stats['winner'].sum() != 1:
                logging.error('No winner {0}/{1}'.format(self.fieldid,
                                                         stats['sourceID'][0]))
                logging.error(stats)
            else:
                logging.debug(stats)

        return stats


###################
# MAIN EXECUTION
###################

if __name__ == "__main__":


    # Which longitude range to process?
    if len(sys.argv) > 1:
        lon1 = int(sys.argv[1])
        lon2 = lon1 + 10
    else:
        raise Exception('Missing longitude strip argument')

    cond_strip = (IPHASQC['is_pdr'] 
                  & (IPHASQC['l'] > lon1-1) 
                  & (IPHASQC['l'] < lon2+1))

    cond_strip = (IPHASQC['id'] == '3680o_dec2003')

    for idx in np.argsort(IPHASQC['seeing_max']):
        if cond_strip[idx]:

            logging.info('Now seaming {0}'.format(IPHASQC['id'][idx]))

            try:
                s = Seamer(IPHASQC['id'][idx], IPHASQC['ra'][idx], IPHASQC['dec'][idx])
                s.run()
            except SeamingException, e:
                logging.error(e)
            
    """
    TABLE = '107.75.fits'
    match = fits.getdata(TABLE, 1)
    for i, row in enumerate(match[-2:-1]):  # range(f[1].size):
        s = Seamer(row)

        if i > 5:
            break
    """
