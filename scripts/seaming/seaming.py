#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Find the primary detections in the IPHAS catalogue."""
from __future__ import division, print_function, unicode_literals
import os
import sys
import numpy as np
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
    DATADIR = "/home/gb/tmp/iphas-dr2/bandmerged"
    # Where to write the output catalogues?
    DESTINATION = "/home/gb/tmp/iphas-dr2/seamed"
else:  # production
    DATADIR = "/car-data/gb/iphas-dr2/bandmerged"
    DESTINATION = "/car-data/gb/iphas-dr2/seamed"

# Where to store temporary files
TMPDIR = '/dev/shm'

# Dir of this script
SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))

# How to execute stilts?
STILTS = 'nice java -Xmx500M -XX:+UseConcMarkSweepGC -jar {0}'.format(
                                 os.path.join(SCRIPTDIR, '../lib/stilts.jar'))

# Where is the IPHAS quality control table?
IPHASQC = fits.getdata('/home/gb/dev/iphas-qc/qcdata/iphas-qc.fits', 1)

# Fields within this radius will be considered to overlap
FIELD_MAXDIST = 0.8  # degrees

# CACHE registers sourceID's for wich a primaryID has already been assigned
CACHE = {}


###########
# CLASSES
###########

class SeamingException(Exception):
    """Exception raised for unexpected problems with a single field."""
    pass


class Seamer(object):
    """Adds the primaryID column to a field's bandmerged catalogue.

    Usage:
       Seamer.run()
    """

    def __init__(self, fieldid, ra, dec, strip):
        """Constructor.

        fieldid -- e.g. 0001_aug2004
        ra -- field R.A. (degrees)
        dec -- field declination (degrees)
        strip -- longitude strip being seamed.
        """
        self.fieldid = fieldid
        self.ra = ra
        self.dec = dec
        self.strip = strip
        # Where to store temporary files?
        self.crossmatch_file = os.path.join(TMPDIR,
                                            'seamtmp_{0}.fits'.format(fieldid))
        self.primaryid_file = os.path.join(TMPDIR,
                                           'seamids_{0}.fits'.format(fieldid))
        # Where to store the results?
        self.output_dir = os.path.join(DESTINATION, 'strip{0}'.format(strip))
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        self.output_file = os.path.join(self.output_dir,
                                        '{0}.fits'.format(fieldid))

    def run(self):
        """Main function."""
        self.overlaps = self.get_overlaps()
        self.crossmatch()
        sourceID, primaryID = self.get_primary_ids()
        self.write_results(sourceID, primaryID)
        self.clean()

    def clean(self):
        """Removes temporary files created by run()."""
        os.remove(self.crossmatch_file)
        os.remove(self.primaryid_file)

    def write_results(self, sourceID, primaryID):
        """Writes a new catalogue with primaryID to disk."""
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
        log.debug(stilts_cmd)
        status = os.system(stilts_cmd)
        log.info('Adding column: '+str(status))
        return status

    def get_overlaps(self):
        """Returns the list of fields which overlap."""
        ddec = self.dec - IPHASQC['dec']
        dra = (self.ra - IPHASQC['ra']) * np.cos(np.radians(self.dec))
        d = (ddec ** 2 + dra ** 2) ** 0.5
        idx = (IPHASQC['is_pdr']
               & (d < FIELD_MAXDIST)
               & (self.fieldid != IPHASQC['id']))
        return IPHASQC['id'][idx]

    def filename(self, fieldid):
        """Returns the path of the bandmerged catalogue for a given field."""
        return os.path.join(DATADIR, '{0}.fits').format(fieldid)

    def get_crossmatch_command(self):
        """Return the stilts command to crossmatch overlapping fields."""
        # Operations to perform on all tables
        icmd = 'addcol bands "sum(array(NULL_r?0:1,NULL_i?0:1,NULL_ha?0:1))"; '
        icmd += 'keepcols "sourceID fieldID ra dec bands errBits seeing rAxis"'
        # Keywords in stilts command
        config = {'STILTS': STILTS,
                  'NIN': len(self.overlaps) + 1,
                  'IN1': self.filename(self.fieldid),
                  'ICMD': icmd,
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
        """Carry out the crossmatching of overlapping fields."""
        cmd = self.get_crossmatch_command()
        log.debug('Overlaps: {0}'.format(self.overlaps))
        log.debug(cmd)
        status = os.system(cmd)
        log.info('Crossmatch: '+str(status))
        if status == 256:
            raise SeamingException('Crossmatch failed for %s' % self.fieldid)
        return status

    def get_primary_ids(self):
        """Returns the tuple (sourceID,primarySourceID)."""
        # Open the file with the sources crossmatched across fields
        try:
            matchtable = fits.getdata(self.crossmatch_file, 1)
        except IOError, e:
            raise SeamingException(e)

        idx_all = (np.arange(self.overlaps.size+1) + 1)  # 1 2 3 ...
        sourceid_keys = ['sourceID_{0}'.format(i) for i in idx_all]

        primary_ids = []  # will hold the key result
        for source in matchtable:
            mySourceID = source['sourceID_1']
            # If we have encountered the source before, we enforce the
            # previous conclusion to ensure consistency
            if mySourceID in CACHE:
                winnerID = CACHE[mySourceID]
            else:

                # Load the sourceIDs of the competitors
                ids = np.array([source[key] for key in sourceid_keys])
                # Filter out the missing IDs (i.e. negative values in FITS)
                idx = idx_all[ids > 0]

                if len(idx) == 1:  # No alternatives
                    winnerID = mySourceID
                    CACHE[mySourceID] = winnerID
                else:

                    my_sourceids = ids[ids > 0]
                    #my_fields = np.array([source['fieldID_%s'%i] for i in idx])
                    my_bands = np.array([source['bands_%s'%i] for i in idx])
                    my_errbits = np.array([source['errBits_%s'%i] for i in idx])
                    my_seeing = np.array([source['seeing_%s'%i] for i in idx])
                    my_raxis = np.array([source['rAxis_%s'%i] for i in idx])

                    winner = self.get_winner(my_bands, my_errbits,
                                             my_seeing, my_raxis)
                    winnerID = my_sourceids[winner]

                    for candidate in my_sourceids:
                        CACHE[candidate] = winnerID
                    #print(stats)
                    #break

            primary_ids.append(winnerID)

        return (matchtable['sourceID_1'], np.array(primary_ids))

    def get_winner(self, bands, errbits, seeing, raxis):
        """
        Returns the index of the primary source given stats.

        The ranking criteria are
          1. prefer filter coverage (3 > 2 > 1);
          2. prefer least severe error flags;
          3. prefer best seeing amongst the filters
             (but tolerate up to 20% of the best value);
          4. prefer smallest distance from optical axis.
        """
        # This array keeps track of which candidates are still in the running
        winning = np.array([True for i in bands])

        # Discard source with less bands
        idx_loser = (bands < np.max(bands))
        winning[idx_loser] = False

        # Discard sources with large errBits
        idx_loser = (errbits > np.min(errbits[winning]))
        winning[idx_loser] = False

        # Discard sources with poor seeing
        idx_loser = (seeing > 1.2*np.min(seeing[winning]))
        winning[idx_loser] = False

        # Finally, the winner is the one closest to the optical axis
        idx_loser = (raxis > np.min(raxis[winning]))
        winning[idx_loser] = False

        # This shouldn't happen
        if winning.sum() != 1:
            log.warning('Object w/o clear winner in {0}'.format(self.fieldid))

        # Return the index of the winner
        return winning.nonzero()[0][0]


###########
# FUNCTIONS
###########

def run_strip(strip):
    # Strips are defined by the start longitude of a 10 deg-wide strip
    assert(strip in np.arange(30, 210+1, 10))
    # So which are our boundaries?
    # Note: we must allow an extree for border overlaps!
    lon1 = strip - 1
    lon2 = lon1 + 10 + 1
    cond_strip = (IPHASQC['is_pdr']
                  & (IPHASQC['l'] > lon1)
                  & (IPHASQC['l'] < lon2))

    #cond_strip = (IPHASQC['id'] == '3693_dec2003')

    # Seam fields; do the best-seeing fields first!
    for idx in np.argsort(IPHASQC['seeing_max']):
        if cond_strip[idx]:

            log.info('Strip {0}: seaming {1}'.format(strip, IPHASQC['id'][idx]))

            try:
                s = Seamer(IPHASQC['id'][idx],
                           IPHASQC['ra'][idx],
                           IPHASQC['dec'][idx],
                           strip)
                s.run()
            except SeamingException, e:
                log.error(str(e))
            except Exception, e:
                log.error('%s: *UNEXPECTED EXCEPTION*: %s' % (IPHASQC['id'][idx], e))


def run_all(lon1=30, lon2=210, ncores=4):
    """ Band-merge all fields """
    log.info('Seaming in longitude strips {0}-{1}'.format(lon1, lon2))
    strips = np.arange(lon1, lon2+1, 10)

    # Distribute the work over ncores
    p = Pool(processes=ncores)
    results = p.imap(run_strip, strips)
    for i in results:
        pass


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
