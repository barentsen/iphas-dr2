#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fix the World Coordinate System (WCS) for IPHAS fields with known problems.

This script takes a list of IPHAS runs which are known to have a poor
astrometric WCS solution. The script then provides an interactive display
to allow the astrometry to be compared against 2MASS, and allows the WCS
parameters to be modified using buttons on the keyboard.
"""
from __future__ import division, print_function
import numpy as np
import os
import urllib2
import tempfile
from astropy.io import fits
from astropy.io import votable
from astropy.io import ascii
from astropy import wcs
from astropy import log
from matplotlib import pyplot as plt

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Hywel Farnhill', 'Robert Greimel', 'Janet Drew']


################################
# CONSTANTS & CONFIGURATION
################################

DATADIR = '/media/0133d764-0bfe-4007-a9cc-a7b1f61c4d1d/iphas'
EXTS = [1, 2, 3, 4]  # Which extensions to expect in the fits catalogues?
PXSCALE = 0.333  # Arcsec/pix of the CCD


############
# FUNCTIONS
############

def twomass_conesearch(ra, dec, radius):
    """
    Returns arrays of ra and dec from a 2MASS conesearch
    """
    url = 'http://vo.astronet.ru/sai_cas/conesearch?cat=twomass&tab=psc'
    url += '&RA={0}&DEC={1}&SR={2}'.format(ra, dec, radius)
    response = urllib2.urlopen(url)
    tf = tempfile.NamedTemporaryFile()
    tf.write(response.read())
    table = votable.parse_single_table(tf)
    return table.array['ra'], table.array['dec']


##########
# CLASSES
##########

class SurveyExposure():
    """INT Wide Field Camera Exposure"""

    def __init__(self, filename, ccd):
        self.path = os.path.join(DATADIR, filename)
        self.ccd = ccd
        self.fits = fits.open(self.path)
        self.sanitize()

    def sanitize(self):
        """Free the FITS header from known errors and inconsistencies"""
        # Early versions of CASU catalogues chave multiple columns 'Blank'
        # Numpy will throw an exception if multiple columns have the same
        # name, so we need to rename these columns.
        n_columns = len(self.fits[self.ccd].columns)
        for col in range(26, n_columns, 1):
            name = self.fits[self.ccd].columns[col].name
            if name == 'Blank':
                self.fits[self.ccd].columns[col].name = 'Blank%d' % col

        # The headers contain a combination of old- and modern-
        # style WCS parameters for the ZPN projection coefficients, which
        # confuses libwcs. Moreover, in a few cases the keyword values
        # are plainly wrong. Hence we remove the keywords.
        for kw in ['PV1_0', 'PV1_1', 'PV1_2', 'PV1_3',
                   'PV2_0', 'PV2_1', 'PV2_2', 'PV2_3',
                   'PV3_0', 'PV3_1', 'PV3_3', 'PV3_3',
                   'PROJP1', 'PROJP3', 'WAT1_001', 'WAT2_001',
                   'RADECSYS']:
            del self.fits[self.ccd].header[kw]

        # ..and enforce the parameters wich have been used by the pipeline
        self.fits[self.ccd].header['EQUINOX'] = 2000.0
        self.fits[self.ccd].header['PV2_1'] = 1.0
        self.fits[self.ccd].header['PV2_3'] = 220.0
        self.fits[self.ccd].header['CUNIT1'] = 'deg'
        self.fits[self.ccd].header['CUNIT2'] = 'deg'
        self.fits[self.ccd].header['RADESYSa'] = 'ICRS'

    def set_wcs(self, ccd_ref):
        for kw in ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                   'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
            self.fits[self.ccd].header[kw] = ccd_ref.fits[self.ccd].header[kw]

    def column(self, name):
        return self.fits[self.ccd].data.field(name)

    def radec(self, shift_crpix1=0.0, shift_crpix2=0.0,
              shift_cd1_1=0.0, shift_cd1_2=0.0,
              shift_cd2_1=0.0, shift_cd2_2=0.0):
        """Returns ra, dec coordinates for the given ccd"""
        header = dict(self.fits[self.ccd].header)
        header['CRPIX1'] += shift_crpix1
        header['CRPIX2'] += shift_crpix2
        header['CD1_1'] += shift_cd1_1
        header['CD1_2'] += shift_cd1_2
        header['CD2_1'] += shift_cd2_1
        header['CD2_2'] += shift_cd2_2
        mywcs = wcs.WCS(header, relax=True)
        # Only return stars and extended objets
        mask = ((self.column('Classification') == -1) 
                | (self.column('Classification') == -2)
                | (self.column('Classification') == 1) 
                | (self.column('Classification') == -3))
        ra, dec = mywcs.wcs_pix2world(self.column('X_coordinate')[mask],
                                      self.column('Y_coordinate')[mask],
                                      1)
        return (ra, dec)


class WCSTuner():
    """
    Class to create a display of the bad vs reference coordinates;
    allowing the WCS to be tuned using keyboard shortcuts.
    """

    def __init__(self, filename_bad, filename_ref, ccd):
        self.ccd = ccd
        self.filename_bad = filename_bad
        self.filename_ref = filename_ref
        self.run_bad = SurveyExposure(filename_bad, ccd)
        self.run_ref = SurveyExposure(filename_ref, ccd)

        self.fig = plt.figure(figsize=(22, 12))
        self.fig.suptitle('Fixing {3} (CCD #{2}) using {4}\n{0} vs {1}'.format(
                          filename_bad, filename_ref,
                          ccd, self.badhdr('OBJECT'), self.refhdr('OBJECT')),
                          fontsize=20)
        self.fig.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.92,
                                 wspace=0.01, hspace=0.01)
        self.fig.canvas.mpl_connect('key_press_event', self.keypress)
        self.axes = [self.fig.add_subplot(321),
                     self.fig.add_subplot(322),
                     self.fig.add_subplot(312),
                     self.fig.add_subplot(325),
                     self.fig.add_subplot(326)]
        for ax in self.axes:
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

        self.shift_crpix1 = 0.0
        self.shift_crpix2 = 0.0
        self.shift_cd1_1, self.shift_cd1_2 = 0.0, 0.0
        self.shift_cd2_1, self.shift_cd2_2 = 0.0, 0.0

    def keypress(self, event):
        """Adjusts WCS parameters by responding to keyboard events.

        arrows: move CRPIX
        []/'  : move CRPIX (fast)
        1/2   : increase/decrease CD1_1
        3/4   : increase/decrease CD1_2
        5/6   : increase/decrease CD2_1
        7/8   : increase/decrease CD2_2
        r     : reset WCS
        w     : write results
        """
        stepsize = 0.15  # pixels, i.e. 0.05 arcsec
        if event.key == 'right':
            if self.ccd == 2:
                self.shift_crpix1 -= stepsize
            else:
                self.shift_crpix2 += stepsize
        elif event.key == 'left':
            if self.ccd == 2:
                self.shift_crpix1 += stepsize
            else:
                self.shift_crpix2 -= stepsize
        elif event.key == ']':
            if self.ccd == 2:
                self.shift_crpix1 -= 5*stepsize
            else:
                self.shift_crpix2 += 5*stepsize
        elif event.key == '[':
            if self.ccd == 2:
                self.shift_crpix1 += 5*stepsize
            else:
                self.shift_crpix2 -= 5*stepsize
        elif event.key == 'up':
            if self.ccd == 2:
                self.shift_crpix2 += stepsize
            else:
                self.shift_crpix1 += stepsize
        elif event.key == 'down':
            if self.ccd == 2:
                self.shift_crpix2 -= stepsize
            else:
                self.shift_crpix1 -= stepsize
        elif event.key == "'":
            if self.ccd == 2:
                self.shift_crpix2 += 5*stepsize
            else:
                self.shift_crpix1 += 5*stepsize
        elif event.key == '/':
            if self.ccd == 2:
                self.shift_crpix2 -= 5*stepsize
            else:
                self.shift_crpix1 -= 5*stepsize
        elif event.key == '1':
            self.shift_cd1_1 += 1e-8
        elif event.key == '2':
            self.shift_cd1_1 -= 1e-8
        elif event.key == '3':
            self.shift_cd1_2 += 1e-8
        elif event.key == '4':
            self.shift_cd1_2 -= 1e-8
        elif event.key == '5':
            self.shift_cd2_1 += 1e-8
        elif event.key == '6':
            self.shift_cd2_1 -= 1e-8
        elif event.key == '7':
            self.shift_cd2_2 += 2e-8
        elif event.key == '8':
            self.shift_cd2_2 -= 2e-8
        elif event.key == 'r':  # Reset
            self.shift_crpix1 = 0.0
            self.shift_crpix2 = 0.0
            self.shift_cd1_1 = 0.0
            self.shift_cd1_2 = 0.0
            self.shift_cd2_1 = 0.0
            self.shift_cd2_2 = 0.0
        elif event.key == 'o':  # OK, field does not need fix
            self.mark_done()
        elif event.key == 'w':  # Write WCS
            self.write()  # Save the results to a csv file
            self.mark_done()
        # Now update the plot for the slight change in WCS parameters
        self.update()

    def refhdr(self, keyword):
        return self.run_ref.fits[self.ccd].header[keyword]

    def badhdr(self, keyword):
        return self.run_bad.fits[self.ccd].header[keyword]

    def mark_done(self):
        filename_done = 'runs-done.tbl'
        done = ''
        if not os.path.exists(filename_done):
            done += 'catalogue\n'  # Header line
        done += self.filename_bad+'\n'
        out = open(filename_done, 'a')
        out.write(done)
        out.close()
        log.debug('Field has been marked done')

    def write(self):
        filename = 'wcs-fixes.csv'
        csv = ''
        if not os.path.exists(filename):
            csv += 'RUN,CCD,CRVAL1,CRVAL2,CRPIX1,CRPIX2,'
            csv += 'CD1_1,CD1_2,CD2_1,CD2_2,REFRUN\n'

        values = []
        values.append(str(self.badhdr('RUN')))
        values.append(str(self.ccd))
        values.append(str(self.refhdr('CRVAL1')))
        values.append(str(self.refhdr('CRVAL2')))
        values.append(str(self.refhdr('CRPIX1') + self.shift_crpix1))
        values.append(str(self.refhdr('CRPIX2') + self.shift_crpix2))
        values.append(str(self.refhdr('CD1_1') + self.shift_cd1_1))
        values.append(str(self.refhdr('CD1_2') + self.shift_cd1_2))
        values.append(str(self.refhdr('CD2_1') + self.shift_cd2_1))
        values.append(str(self.refhdr('CD2_2') + self.shift_cd2_2))
        values.append(str(self.refhdr('RUN')))
        csv += ','.join(values)

        out = open(filename, 'a')
        out.write(csv+'\n')
        out.close()

        log.info("Wrote to %s:\n%s" % (filename, csv))

    def mask(self, ra, dec, axnum):
        idx = ((ra > self.xlim[axnum][0]) & (ra < self.xlim[axnum][1])
               & (dec > self.ylim[axnum][0]) & (dec < self.ylim[axnum][1]))
        return ra[idx], dec[idx]

    def start(self):
        log.debug('Preparing to show %s' % self.filename_bad)
        ra, dec = self.run_ref.radec()
        bad_ra, bad_dec = self.run_bad.radec()
        tm_ra, tm_dec = twomass_conesearch(np.median(ra), np.median(dec), 0.25)

        log.debug('Found %s 2MASS sources' % len(tm_ra))

        lim = 60/3600.0
        self.xlim = [[ra.min(), ra.min()+2*lim],
                     [ra.max()-2*lim, ra.max()],
                     [np.median(ra)-2*lim, np.median(ra)+2*lim],
                     [ra.min(), ra.min()+2*lim],
                     [ra.max()-2*lim, ra.max()]]

        self.ylim = [[dec.max()-2*lim, dec.max()],
                     [dec.max()-2*lim, dec.max()],
                     [np.median(dec)-1.5*lim, np.median(dec)+1.5*lim],
                     [dec.min(), dec.min()+2*lim],
                     [dec.min(), dec.min()+2*lim]]

        for i, ax in enumerate(self.axes):
            myra, mydec = self.mask(bad_ra, bad_dec, i)
            ax.scatter(myra, mydec, marker='o', facecolor='None',
                       edgecolor='#888888', linewidth=0.5, s=10)

            # Plot reference CCD as red circles
            myra, mydec = self.mask(ra, dec, i)
            ax.scatter(myra, mydec, marker='o', facecolor='None',
                       edgecolor='red', linewidth=0.5, s=10)

            # Plot 2MASS as blue crosses
            myra, mydec = self.mask(tm_ra, tm_dec, i)
            ax.scatter(tm_ra, tm_dec, marker='x', facecolor='white',
                       edgecolor='blue', linewidth=0.5, s=90)

            ax.set_xlim(self.xlim[i])
            ax.set_ylim(self.ylim[i])

        self.run_bad.set_wcs(self.run_ref)

        self.update()
        plt.show()

    def update(self):
        if hasattr(self, 'crosses'):
            for crosses in self.crosses:
                crosses.set_visible(False)
                del crosses

        self.crosses = []

        ra, dec = self.run_bad.radec(self.shift_crpix1, self.shift_crpix2,
                                     self.shift_cd1_1, self.shift_cd1_2,
                                     self.shift_cd2_1, self.shift_cd2_2)
        for i, ax in enumerate(self.axes):
            myra, mydec = self.mask(ra, dec, i)
            mycrosses = ax.scatter(myra, mydec, marker='+', facecolor=None,
                                   edgecolor='black', linewidth=0.5, s=200)
            self.crosses.append(mycrosses)

        plt.draw()


#######
# MAIN
#######

if __name__ == '__main__':
    log.setLevel('DEBUG')

    done = ascii.read('wcs-fixes.csv')  # File listing runs already fixed
    #for row in ascii.read('runs-todo.tbl'):
    for row in ascii.read('needfix.txt'):
        # Have we already corrected this run/ccd combination?
        myrun = int(row['catalogue'].split('/')[1].split('_')[0][1:])
        if (myrun in done['RUN']
                and row['ccd'] in done['CCD'][done['RUN'] == myrun]):
            log.info('Skipping {0} CCD#{1}: already done'.format(myrun,
                                                                 row['ccd']))
            continue

        tuner = WCSTuner(row['catalogue'], row['ref'], row['ccd'])
        tuner.start()
    log.info('All is done.')
