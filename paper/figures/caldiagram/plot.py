"""Plots the calibrated and uncalibrated CCD over a large area."""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import log
from dr2 import constants

fields = constants.IPHASQC[constants.IPHASQC_COND_RELEASE]
mask = (fields['l'] > 160.0) & (fields['l'] < 200.)
ids = fields['id'][mask]
log.info('Plotting {0} fields.'.format(len(ids)))

calibrated = True

if calibrated:
    CATALOGUE_PATH = '/home/gb/tmp/iphas-dr2-rc6/bandmerged-calibrated/'
else:
    CATALOGUE_PATH = '/home/gb/tmp/iphas-dr2-rc6/bandmerged/'

COLORMAP = matplotlib.colors.LinearSegmentedColormap.from_list('mymap', 
    ['#bd0026', '#f03b20', '#fd8d3c', '#fecc5c', '#ffffb2'])



rmi, rmha = [], []

for field in ids:
    filename = os.path.join(CATALOGUE_PATH, field+'.fits')
    log.info(filename)
    d = fits.getdata(filename, 1)
    veryreliable = (d['reliable']
                    & (d['pStar'] > 0.9)
                    & -d['deblend']
                    & -d['brightNeighb']
                    & (d['r'] < 18)
                    )
    rmi.append(d['rmi'][veryreliable])
    rmha.append(d['rmha'][veryreliable])
    #plt.scatter(d['rmi'][veryreliable], d['rmha'][veryreliable])


plt.figure()
plt.subplots_adjust(0.16, 0.20, 0.98, 0.98)

vmax = 500
hist, xedges, yedges = np.histogram2d(np.concatenate(rmi),
                                      np.concatenate(rmha),
                     range=[[-0.3, 3.1], [-0.1, 2.1]], 
                      bins=[320, 300])
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
hist[hist == 0] = np.nan
plt.imshow(hist.T, extent=extent, 
           interpolation='nearest',
           origin='lower',
           aspect='auto',
           vmin=0,
           vmax=vmax,
           cmap=COLORMAP)

if calibrated:
    label = '(b) After re-calibration'
else:
    label = '(a) Before re-calibration'
plt.text(0.06, 0.92, label,
         horizontalalignment='left',
         verticalalignment='top',
         transform=plt.gca().axes.transAxes,
         fontsize=11)

plt.xlabel('$r$ - $i$')
plt.ylabel('$r$ - $\\rm H\\alpha$')
plt.xlim([-0.3, 3.])
plt.ylim([-0.1, 1.5])

if calibrated:
    output = 'ccd-calibrated.pdf'
else:
    output = 'ccd-uncalibrated.pdf'
plt.savefig(output, dpi=300)
plt.close()
