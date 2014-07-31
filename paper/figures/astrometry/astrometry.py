"""Plots a histogram showing positional residuals between IPHAS and UCAC4."""
import pylab as p
from astropy.table import Table

t = Table.read('iphas-x-ucac4.fits')

p.figure()
p.hist(t['Separation']*3600., range=[0.0, 0.5], bins=50, lw=1.0, color='#dddddd')
p.minorticks_on()
p.xlim([0,0.5])
p.xlabel('Astrometric residuals against UCAC4 [arcsec]', fontsize=9)
p.ylabel('Stars', fontsize=9)
p.tight_layout(pad=0.5)
p.savefig('astrometry.pdf')
