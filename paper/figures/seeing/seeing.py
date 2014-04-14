from astropy.io import fits
import pylab as p

d = fits.getdata('/home/gb/dev/iphas-qc/qcdata/iphas-qc.fits', 1)

mask = d['is_dr2']


for band in ['r', 'i', 'ha']:

    p.figure()
    p.minorticks_on()

    n, bins, patches = p.hist(d['seeing_'+band][mask],
                              range=[0.5, 2.5],
                              bins=20, 
                              #histtype='stepfilled',
                              #cumulative=True,
                              lw=1.0,
                              color='#dddddd')

    p.ylabel('Fields')
    p.xlabel("PSF FWHM [arcsec]")
    p.xlim([0.5, 2.5])

    if band == 'ha':
        label = '$\\rm H\\alpha$'
    else:
        label = '$'+band+'$'

    p.text(0.87, 0.9, label,
        horizontalalignment='center',
        verticalalignment='top',
        transform=p.axes().transAxes,
        fontsize=16)

    p.tight_layout(pad=0.5)
    p.savefig('seeing_{0}.pdf'.format(band))
    p.ylim(0, 2500)

