"""Plots example colour/magnitude diagrams along different sightlines."""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib import colors
import numpy as np
from astropy.io import ascii
from astropy import log

from dr2.sql import sql
db = sql.SourceCatalogueDB()

COLORMAP = colors.LinearSegmentedColormap.from_list('mymap', 
    ['#bd0026', '#f03b20', '#fd8d3c', '#fecc5c', '#ffffb2'])

# Tracks from Drew et al. 2005
drew2005_a0 = {'rmi': [0.029, 0.699, 1.352, 1.991, 2.616],
               'rmha': [0.001, 0.199, 0.355, 0.468, 0.544]}

drew2005_dwarfs = {'ebv0': {
                    'rmi': [-0.176, 0.029, 0.212,
                            0.445, 0.650, 0.903, 1.144,
                            1.521], #, 1.829], #, 2.514],
                    'rmha': [0.067, 0.001, 0.114,
                             0.278, 0.390, 0.499, 0.624,
                             0.820]} #, 0.889]}#, 1.060]}
                    }

drew2005_giants = {'ebv0': {
                    'rmi': [-0.180, 0.013, 0.177, 0.378, 
                            0.505, 0.808, 0.935, 1.040,
                            1.311, 1.652], #, 1.941],
                    'rmha': [0.075, 0.022, 0.103, 0.236,
                             0.284, 0.381, 0.473, 0.512,
                             0.564, 0.644]}, #, 0.725]},
                   'ebv1': {
                    'rmi': [0.489, 0.682, 0.844, 1.035,
                            1.157, 1.456, 1.581, 1.679,
                            1.957, 2.296, 2.581],
                    'rmha': [0.282, 0.220, 0.290, 0.408,
                             0.445, 0.525, 0.611, 0.641,
                             0.687, 0.753, 0.821]},
                   'ebv2': {
                    'rmi': [1.142, 1.335, 1.497, 1.677,
                            1.793, 2.088, 2.212, 2.304,
                            2.586, 2.924, 3.205],
                    'rmha': [0.446, 0.375, 0.435, 0.537,
                             0.564, 0.629, 0.708, 0.730,
                             0.770, 0.822, 0.879]}
                }


def plot_ccd(l, b=0):
    data = db.conesearch_gal(l, b, 
                             radius=np.sqrt(1/np.pi)*3600,
                             condition='AND reliable '
                                       'AND NOT deblend '
                                       'AND NOT brightNeighb '
                                       'AND pStar > 0.9 '
                                       'AND rErr < 0.05 '
                                       'AND iErr < 0.05 '
                                       'AND haErr < 0.05')
    log.info('Found {0} sourcesat l={1}'.format(data.size, l))    

    plt.figure()
    ax = plt.subplot()
    plt.subplots_adjust(0.16, 0.20, 0.98, 0.98)

    vmax = 25
    hist, xedges, yedges = np.histogram2d(data['rmi'],
                                          data['rmha'], 
                         range=[[-0.1, 3.1], [-0.1, 2.1]], 
                          bins=[300, 300])
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
    hist[hist == 0] = np.nan
    plt.imshow(hist.T, extent=extent, 
               interpolation='nearest',
               origin='lower',
               aspect='auto',
               vmin=0,
               vmax=vmax,
               cmap=COLORMAP)

    minorLocator = MultipleLocator(0.1)
    ax.xaxis.set_minor_locator(minorLocator)
    minorLocator2 = MultipleLocator(0.05)
    ax.yaxis.set_minor_locator(minorLocator2)

    if l == 180:
        plt.annotate('A0V', (0.029, 0.00),
                     xytext=(-5, 15), textcoords='offset points', ha="center", fontsize=6,
                     arrowprops=dict(arrowstyle="-", connectionstyle="arc", shrinkA=0, shrinkB=0))

        plt.annotate('K0V', (0.445, 0.278),
                     xytext=(-15, 10), textcoords='offset points', ha="center", fontsize=6,
                     arrowprops=dict(arrowstyle="-", connectionstyle="arc", shrinkA=0, shrinkB=0))

        plt.annotate('M0V', (0.903, 0.499),
                     xytext=(-15, 10), textcoords='offset points', ha="center", fontsize=6,
                     arrowprops=dict(arrowstyle="-", connectionstyle="arc", shrinkA=0, shrinkB=0))

        plt.annotate('M3V', (1.521, 0.820),
                     xytext=(-15, 10), textcoords='offset points', ha="center", fontsize=6,
                     arrowprops=dict(arrowstyle="-", connectionstyle="arc", shrinkA=0, shrinkB=0))

        plt.annotate('M4III', (1.652, 0.644),
                     xytext=(15, -10), textcoords='offset points', ha="center", fontsize=6,
                     arrowprops=dict(arrowstyle="-", connectionstyle="arc", shrinkA=0, shrinkB=0))

        plt.annotate('E(B-V)=1', (0.699, 0.199),
                     xytext=(10, -12), textcoords='offset points', ha="center", fontsize=6,
                     arrowprops=dict(arrowstyle="-", connectionstyle="arc", shrinkA=0, shrinkB=0))
        plt.annotate('E(B-V)=2', (1.352, 0.355),
                     xytext=(10, -12), textcoords='offset points', ha="center", fontsize=6,
                     arrowprops=dict(arrowstyle="-", connectionstyle="arc", shrinkA=0, shrinkB=0))
        plt.annotate('E(B-V)=3', (1.991, 0.468),
                     xytext=(10, -12), textcoords='offset points', ha="center", fontsize=6,
                     arrowprops=dict(arrowstyle="-", connectionstyle="arc", shrinkA=0, shrinkB=0))

    plt.plot(drew2005_a0['rmi'][0:4],
             drew2005_a0['rmha'][0:4],
             lw=2.0, ls='dashed', c='black')
    plt.plot(drew2005_dwarfs['ebv0']['rmi'],
             drew2005_dwarfs['ebv0']['rmha'],
             lw=1.5, c='black')
    plt.plot(drew2005_giants['ebv0']['rmi'],
             drew2005_giants['ebv0']['rmha'],
             lw=2.5, c='black')

    plt.xlabel('$r$ - $i$')
    plt.ylabel('$r$ - $\\rm H\\alpha$')
    plt.xlim([-0.3, 3.])
    plt.ylim([-0.1, 1.3])

    plt.text(0.06, 0.92, '$(\ell,b)=({0}^\circ, {1}^\circ)$'.format(l, b),
             horizontalalignment='left',
             verticalalignment='top',
             transform=ax.transAxes,
             fontsize=11)

    plt.savefig('ccd-{0}-{1}.pdf'.format(l, b), dpi=200)
    plt.close()


def plot_cmd(l, b=0):
    data = db.conesearch_gal(l, b, 
                             radius=np.sqrt(1/np.pi)*3600,
                             condition='AND reliable '
                                       'AND NOT deblend '
                                       'AND NOT brightNeighb '
                                       'AND pStar > 0.9 '
                                       'AND rErr < 0.05 '
                                       'AND iErr < 0.05 '
                                       'AND haErr < 0.05'
                                       )
    log.info('Found {0} sourcesat l={1}'.format(data.size, l))    

    plt.figure()
    ax = plt.subplot()
    plt.subplots_adjust(0.16, 0.20, 0.98, 0.98)

    vmax = 25
    hist, xedges, yedges = np.histogram2d(data['rmi'],
                                          data['r'], 
                         range=[[-0.1, 3.1], [12, 23]], 
                          bins=[300, 300])
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
    hist[hist == 0] = np.nan
    plt.imshow(hist.T, extent=extent, 
               interpolation='nearest',
               origin='lower',
               aspect='auto',
               vmin=0,
               vmax=vmax,
               cmap=COLORMAP)

    minorLocator = MultipleLocator(0.1)
    ax.xaxis.set_minor_locator(minorLocator)
    minorLocator2 = MultipleLocator(0.5)
    ax.yaxis.set_minor_locator(minorLocator2)

    if l == 180:
        plt.annotate('E(B-V)=1', (drew2005_a0['rmi'][1], 12 + 0.843*3.1),
                     xytext=(20, 20), textcoords='offset points', ha="center", fontsize=6,
                     arrowprops=dict(arrowstyle="-", connectionstyle="arc", shrinkA=0, shrinkB=0))
        plt.annotate('E(B-V)=2', (drew2005_a0['rmi'][2], 12 + 0.843*3.1*2),
                     xytext=(20, 20), textcoords='offset points', ha="center", fontsize=6,
                     arrowprops=dict(arrowstyle="-", connectionstyle="arc", shrinkA=0, shrinkB=0))
        plt.annotate('E(B-V)=3', (drew2005_a0['rmi'][3], 12 + 0.843*3.1*3),
                     xytext=(20, 20), textcoords='offset points', ha="center", fontsize=6,
                     arrowprops=dict(arrowstyle="-", connectionstyle="arc", shrinkA=0, shrinkB=0))

    plt.plot(drew2005_a0['rmi'][0:5],
             12 + 0.843*3.1*np.arange(len(drew2005_a0['rmi'][0:5])),
             lw=2.0, ls='dashed', c='black')   

    padova = ascii.read('padova-1gyr.dat', header_start=12, data_start=0)
    plt.plot(padova['gR']-padova['gI'],
             padova['gR']+12,
             c='black', lw=1.5)

    plt.xlabel('$r$ - $i$')
    plt.ylabel('$r$')
    plt.xlim([-0.3, 3.])
    plt.ylim([20., 11.5])

    plt.savefig('cmd-{0}-{1}.pdf'.format(l, b), dpi=200)
    plt.close()


if __name__ == '__main__':
    plot_ccd(30, 0)
    plot_ccd(45, 2)
    plot_ccd(180, 3)
    plot_cmd(30, 0)
    plot_cmd(45, 2)
    plot_cmd(180, 3)
