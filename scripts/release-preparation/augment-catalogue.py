"""Adds a few finishing touches to the IPHAS DR2 binary FITS catalogues.

It will add TDISP keywords, sanitise TUNIT keywords, add checksums,
and add origin information.
"""
import numpy as np
from astropy.io import fits
from astropy import log


def augment(filename_origin, filename_target):
    log.info('Opening {0}'.format(filename_origin))
    f = fits.open(filename_origin)

    for i in np.arange(1, f[1].header['TFIELDS']+1, 1):  # Loop over all columns

        # Set an appropriate TDISP keyword for floating points
        name = f[1].header['TTYPE{0}'.format(i)]
        if name in ['ra', 'dec', 'l', 'b']:
            f[1].header['TDISP{0}'.format(i)] = 'F9.5'
        if name in ['posErr', 'pStar', 'pGalaxy', 'pNoise', 'pSaturated',
                    'rGauSig', 'rEll', 'rSeeing',
                    'iGauSig', 'iEll', 'iSeeing',
                    'haGauSig', 'haEll', 'haSeeing',
                    'seeing']:
            f[1].header['TDISP{0}'.format(i)] = 'F4.2'
        if name in ['mergedClassStat', 'rmi', 'rmha',
                    'r', 'rErr', 'rPeakMag', 'rPeakMagErr',
                    'rAperMag1', 'rAperMag1Err', 'rAperMag3', 'rAperMag3Err',
                    'rClassStat',
                    'i', 'iErr', 'iPeakMag', 'iPeakMagErr',
                    'iAperMag1', 'iAperMag1Err', 'iAperMag3', 'iAperMag3Err',
                    'iClassStat', 'iXi', 'iEta',
                    'ha', 'haErr', 'haPeakMag', 'haPeakMagErr',
                    'haAperMag1', 'haAperMag1Err', 'haAperMag3', 'haAperMag3Err',
                    'haClassStat', 'haXi', 'haEta',
                    'r2', 'rErr2', 'i2', 'iErr2', 'ha2', 'haErr2']:
            f[1].header['TDISP{0}'.format(i)] = 'F5.2'
        if name in ['rPa', 'iPa', 'haPa']:
            f[1].header['TDISP{0}'.format(i)] = 'F5.1'
        if name in ['rMJD', 'iMJD', 'haMJD']:
            f[1].header['TDISP{0}'.format(i)] = 'F11.5'
        if name in ['rX', 'rY', 'iX', 'iY', 'haX', 'haY']:
            f[1].header['TDISP{0}'.format(i)] = 'F7.2'

        # Bring unit definitions in line with the FITS standard
        try:
            unit = f[1].header['TUNIT{0}'.format(i)]
            if unit == 'degrees':
                f[1].header['TUNIT{0}'.format(i)] = 'deg'
            if unit == 'Magnitude':
                f[1].header['TUNIT{0}'.format(i)] = 'mag'
            if unit == 'Pixels':
                f[1].header['TUNIT{0}'.format(i)] = 'pixel'
            if unit == 'Arcsec':
                f[1].header['TUNIT{0}'.format(i)] = 'arcsec'

            if unit in ['Sigma', 'Number', 'Flag', 'N-sigma',
                        'bitmask', 'Julian days', 'String']:
                del f[1].header['TUNIT{0}'.format(i)]
        except KeyError:
            pass

    # Make the header more informative
    f[1].header['EXTNAME'] = 'CATALOG'
    f[1].header['ORIGIN'] = 'IPHAS'
    f[1].header['PHOTSYS'] = 'VEGA'
    f[1].header['REFERENC'] = 'Barentsen et al (2014)'
    f[1].header['PRODCATG'] = 'SCIENCE.CATALOGTILE'
    f[1].header['COMMENT'] = ' _____ _____  _    _           _____ '     
    f[1].header['COMMENT'] = '|_   _|  __ \| |  | |   /\    / ____|'
    f[1].header['COMMENT'] = '  | | | |__) | |__| |  /  \  | (___  '
    f[1].header['COMMENT'] = '  | | |  ___/|  __  | / /\ \  \___ \ '
    f[1].header['COMMENT'] = ' _| |_| |    | |  | |/ ____ \ ____) |'
    f[1].header['COMMENT'] = '|_____|_|    |_|  |_/_/    \_\_____/ '
    f[1].header['COMMENT'] = ''
    f[1].header['COMMENT'] = 'This catalogue is part of IPHAS DR2.'
    f[1].header['COMMENT'] = 'For more information, visit http://www.iphas.org.'

    f.writeto(filename_target, checksum=True,
              clobber=True)


if __name__ == '__main__':
    DR2 = '/car-data/gb/iphas-dr2-rc6/concatenated'

    for l in np.arange(25, 220, 5):
        for part in ['a', 'b']:
            origin = DR2+'/light/iphas-dr2-{0:03d}{1}-light.fits'.format(l, part)
            target = DR2+'/light-augmented/iphas-dr2-{0:03d}{1}-light.fits.gz'.format(l, part)
            #augment(origin, target)
            origin = DR2+'/full/iphas-dr2-{0:03d}{1}.fits'.format(l, part)
            target = DR2+'/full-augmented/iphas-dr2-{0:03d}{1}.fits'.format(l, part)
            augment(origin, target)

