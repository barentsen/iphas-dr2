#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Find the primary detections in the IPHAS catalogue."""
from __future__ import division, print_function, unicode_literals
import os
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

FIELD_MAXDIST = 0.8  # degrees

CACHE = {}


###########
# CLASSES
###########

class Seamer(object):

    def __init__(self, fieldid, ra, dec):
        self.fieldid = fieldid
        self.crossmatch_file = '/tmp/seaming_{0}.fits'.format(fieldid)
        self.ra = ra
        self.dec = dec

    def run(self):
        self.overlaps = self.get_overlaps()
        #self.crossmatch()
        primaryID = self.get_primary_ids()
        #print(primaryID)
        self.add_column(primaryID)

    def add_column(self, array):
        col_primaryID = fits.Column(name='primaryID', format='K',
                                    unit='Number', array=array)

        #f = fits.open(self.filename(self.fieldid))
        #f[1].columns.add_col(col_primaryID)
        #f.writeto('/tmp/test2.fits', clobber=True)
        #fits.writeto('/tmp/test2.fits', f[1].data, f[1].header, clobber=True)

        output_filename = '/tmp/test2.fits'

        # Open old file
        f = fits.open(self.filename(self.fieldid), memmap=False)

        tmp = f[1].data['ra']  # deal with crazy bug
        #logging.error(len(array))

        cdefs = f[1].get_coldefs()
        cdefs.add_col(col_primaryID)
        outtabhdu = fits.new_table(cdefs)
        outtabhdu.writeto(output_filename, clobber=True)
        
        #hdu_table = fits.new_table(cols, tbtype='BinTableHDU')
        #fits.new_table(f[1].columns, tbtype='BinTableHDU')

        #hdu_primary = fits.PrimaryHDU()
        #hdulist = fits.HDUList([hdu_primary, hdu_table])
        #f[1].columns.add_col(col_primaryID)
        #f.writeto(output_filename, clobber=True)
        


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
        url = ('http://stri-cluster.herts.ac.uk/~gb/'
               + 'iphas-dr2/iphasSource/{0}.fits')
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
        cmd = "{STILTS} tmatchn matcher=sky params=0.5 "
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
        logging.info('Overlaps: {0}'.format(self.overlaps))
        logging.info(cmd)
        status = os.system(cmd)
        logging.info(status)
        return status

    def get_primary_ids(self):
        d = fits.getdata(self.crossmatch_file, 1)

        primary_ids = []
        for source in d:

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

        return np.array(primary_ids)

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

    for idx in np.argsort(IPHASQC['seeing_max']):
        logging.info('Now seaming {0}'.format(IPHASQC['id'][idx]))

        s = Seamer(IPHASQC['id'][idx], IPHASQC['ra'][idx], IPHASQC['dec'][idx])
        s.run()
        break
    """
    TABLE = '107.75.fits'
    match = fits.getdata(TABLE, 1)
    for i, row in enumerate(match[-2:-1]):  # range(f[1].size):
        s = Seamer(row)

        if i > 5:
            break
    """
