import numpy as np
from astropy.table import Table
from astropy.table import Column
from dr2.constants import IPHASQC


t = Table.read('/home/gb/tmp/dl/iphas-images.fits')

# Run 376022 on the disk received from CASU is a corrupt file
t.remove_row(np.argwhere(t['run'] == 376022)[0][0])
# Run 367744 appeared twice in iphas-qc.fits
t.remove_rows(np.argwhere(t['run'] == 367744)[4:])

# Add the full URL of the image location
urldata = ['www.iphas.org/data/images/'+name[0:4]+'/'+name for name in t['filename']]
url = Column(name='url', data=urldata)
t.add_column(url, 0)
t.remove_column('filename')

# Load auxillary data from the IPHAS-QC file
runs = np.concatenate((IPHASQC['run_r'], IPHASQC['run_i'], IPHASQC['run_ha']))
fields = np.concatenate((IPHASQC['id'], IPHASQC['id'], IPHASQC['id']))
qflags = np.concatenate((IPHASQC['qflag'], IPHASQC['qflag'], IPHASQC['qflag']))
field_dict = dict(zip(runs, fields))
qflag_dict = dict(zip(runs, qflags))

# Add the IPHAS field number
field = Column(name='fieldid', data=[field_dict[r] for r in t['run']])
t.add_column(field)

# Add the DR2 quality grade
qcgrade = Column(name='qcgrade', data=[qflag_dict[r] for r in t['run']])
t.add_column(qcgrade, 3)

t.remove_column('airmass')

t.sort(['run', 'ccd'])

columns = ['run', 'ccd', 'url', 'ra', 'dec', 'band', 'utstart',
           'fieldid', 'qcgrade', 'in_dr2',
           'exptime', 'seeing', 'elliptic', 'skylevel', 'skynoise',
           'photzp', 'confmap',
           'ra_min', 'ra_max', 'dec_min', 'dec_max']
t[columns].write('iphas-images.fits', overwrite=True)

t[columns[2:10]].write('iphas-images.txt', format='ascii.fixed_width')