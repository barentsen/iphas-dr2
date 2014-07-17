"""Plots magnitude distribution."""
import matplotlib.pyplot as plt
from astropy.table import Table

for band in ['r']:
    r = Table.read(band+'.csv', format='ascii.csv')
    r2 = Table.read(band+'_reliable.csv', format='ascii.csv')
    r3 = Table.read(band+'_veryreliable.csv', format='ascii.csv')


    plt.figure(figsize=(3.5, 2.5))

    plt.bar(r[band]-0.25, r['count'], width=0.5, 
            facecolor='#cccccc',
            label='all sources')
    plt.bar(r2[band]-0.25, r2['count'], width=0.5,
            facecolor='#888888',
            label='\emph{a10}')
    plt.bar(r3[band]-0.25, r3['count'], width=0.5,
            facecolor='#444444',
            label='\emph{a10point}')

    plt.legend(loc='upper left', fontsize=9)
    plt.xlim([12.25, 22.75])
    plt.xlabel(band+' magnitude')
    plt.ylabel('Sources')

    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1.))
    #plt.gca().yaxis.set_major_locator(y_majorLocator)

    plt.subplots_adjust(0.15, 0.18, 0.98, 0.9)

    plt.savefig('magdist-'+band+'.pdf')
    plt.close()