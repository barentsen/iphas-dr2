import pylab as p
from astropy.table import Table


t = Table.read('smallmatch.fits')

p.figure()
p.hist(t['angDist'], range=[0.0, 0.5],
                   bins=50, 
                   lw=1.0,
                   color='#dddddd')
p.minorticks_on()
                             

p.xlim([0,0.5])
p.xlabel('{\sc iphas} - {\sc ucac4} position offset [arcsec]')
p.ylabel('Stars')
#tight_layout(pad=0.5)
p.savefig('astrometry.pdf')
