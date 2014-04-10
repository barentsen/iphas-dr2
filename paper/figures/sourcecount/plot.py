"""Plots IPHAS DR2 source counts as a function of longitude."""
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.table import Table

t1 = Table.read('count-all.csv', format='ascii.csv')
t2 = Table.read('count-not-deblend.csv', format='ascii.csv')

mask = (t1['l'] >= 30) & (t1['l'] < 215)

plt.figure(figsize=(7,2.5))
plt.subplots_adjust(0.08, 0.17, 0.98, 0.92)

plt.step(t1['l'][mask]+0.5, t1['count'][mask], 
         linewidth=1.5, color='#2980b9',
         label='All sources')
plt.step(t2['l'][mask]+0.5, t2['count'][mask],
         linewidth=1.0, color='#c0392b',
         label='Non-blended sources')

minorLocator = MultipleLocator(5.)
plt.gca().xaxis.set_minor_locator(minorLocator)
majorLocator = MultipleLocator(20.)
plt.gca().xaxis.set_major_locator(majorLocator)
#minorLocator = MultipleLocator(2.5)
#plt.gca().yaxis.set_minor_locator(minorLocator)

plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))

plt.legend(loc='upper left', fontsize=9)
plt.xlabel('Galactic longitude ($l$)')
plt.ylabel('Source count')
plt.xlim([220, 20])
plt.savefig('sourcecount.pdf')
plt.close()
