from astropy.io import fits
import pylab as p

d = fits.getdata('/home/gb/dev/iphas-qc/qcdata/iphas-qc.fits', 1)

mask = d['is_dr2']


for band in ['r', 'i', 'h']:

    p.figure()
    n, bins, patches = p.hist(d[band+'5sig_judged'][mask], 
                              range=[18,23],
                              bins=50,
                              lw=1.0,
                              color='#dddddd')


    p.minorticks_on()
    p.ylabel('Fields')
    p.xlabel('$5\sigma$ limiting magnitude')
    p.xlim([19, 22.0])
    #p.ylim([0, 14500])
    #p.legend(loc='lower right')

    if band == 'h':
        label = '$\\rm H\\alpha$'
    else:
        label = '$'+band+'$'

    p.text(0.12, 0.9, label,
        horizontalalignment='center',
        verticalalignment='top',
        transform=p.axes().transAxes,
        fontsize=16)

    p.tight_layout(pad=0.5)
    p.ylim(0, 2600)
    p.savefig('depth_{0}.pdf'.format(band))
    

