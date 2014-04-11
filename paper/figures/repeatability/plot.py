"""Plots the mean Poissonion error as a function of magnitude."""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.table import Table

reliable = True

if reliable:
    version = "_errbits0_veryreliable.csv"
else:
    version = "_errbits16.csv"

data_r = Table.read('data/rmr2'+version)
data_i = Table.read('data/imi2'+version)
data_ha = Table.read('data/hamha2'+version)

xlim = [11.25, 22.75]
ylim = [-0.04, 0.49]

majorLocator = MultipleLocator(1.0)
y_majorLocator = MultipleLocator(0.1)

plt.figure(figsize=(3.5, 3.5))
plt.subplots_adjust(0.17, 0.15, 0.95, 0.93, hspace=0, wspace=0)

ax = plt.subplot(3, 1, 1)
plt.plot([10,30], [0.1,0.1], linestyle='dashed', color='#666666', linewidth=0.5, zorder=-5)
plt.plot([10,30], [0.2,0.2], color='#666666', linewidth=0.5, zorder=-5)
plt.errorbar(data_r['mag'], data_r['mean'], yerr=data_r['std'], linestyle='None', marker='o',
             markersize=3.0, c='black', capsize=1.5, elinewidth=1.5)
plt.xlim(xlim)
plt.ylim(ylim)
plt.ylabel('$|r-r_2|$')
plt.gca().xaxis.set_ticklabels([])
plt.gca().xaxis.set_major_locator(majorLocator)
plt.gca().yaxis.set_major_locator(y_majorLocator)

if reliable:
    plt.text(0., 1.05, "(b) \emph{veryReliable} sources", fontsize=11,
             ha='left', va='bottom', transform=ax.transAxes)
else:
    plt.text(0., 1.05, "(a) all sources", fontsize=11,
             ha='left', va='bottom', transform=ax.transAxes)

plt.subplot(3, 1, 2)
plt.plot([10,30], [0.1,0.1], linestyle='dashed', color='#666666', linewidth=0.5)
plt.plot([10,30], [0.2,0.2], color='#666666', linewidth=0.5)
plt.errorbar(data_i['mag'], data_i['mean'], yerr=data_i['std'], linestyle='None', marker='o',
             markersize=3.0, c='black', capsize=1.5, elinewidth=1.5)
plt.xlim(xlim)
plt.ylim(ylim)
plt.ylabel('$|i-i_2|$')
plt.gca().xaxis.set_ticklabels([])
plt.gca().xaxis.set_major_locator(majorLocator)
plt.gca().yaxis.set_major_locator(y_majorLocator)

plt.subplot(3, 1, 3)
plt.plot([10,30], [0.1,0.1], linestyle='dashed', color='#666666', linewidth=0.5)
plt.plot([10,30], [0.2,0.2], color='#666666', linewidth=0.5)
plt.errorbar(data_ha['mag'], data_ha['mean'], yerr=data_ha['std'], linestyle='None', marker='o',
             markersize=3.0, c='black', capsize=1.5, elinewidth=1.5)
plt.xlim(xlim)
plt.ylim(ylim)
plt.ylabel('$|\\mathrm{H\\alpha}-\\mathrm{H\\alpha_2}|$')
plt.xticks([])
plt.gca().xaxis.set_major_locator(majorLocator)
plt.gca().yaxis.set_major_locator(y_majorLocator)

plt.xlabel('Magnitude')

if reliable:
    plt.savefig('pairmag-reliable.pdf')
else:
    plt.savefig('pairmag.pdf')
    
plt.close()
