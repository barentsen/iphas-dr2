"""Script to create a user-friendly index of IPHAS image meta data.
"""
import numpy as np
from astropy.table import Table
from astropy.table import Column
from dr2.constants import IPHASQC

# Index of images found by the DR2 pipeline
# ie. produced by dr2.images.prepare_images()
t = Table.read('iphas-images-pipeline.fits')

# Run 376022 on the disk received from CASU is a corrupt file
t.remove_row(np.argwhere(t['run'] == 376022)[0][0])
# Run 367744 appeared twice in iphas-qc.fits
t.remove_rows(np.argwhere(t['run'] == 367744)[4:])

# Add the URL of the image location
urldata = ['http://www.iphas.org/data/images/'+name[0:4]+'/'+name for name in t['filename']]
url = Column(name='url', data=urldata)
t.add_column(url, 0)
t.remove_column('filename')

# Load auxillary data from the IPHAS-QC file
runs = np.concatenate((IPHASQC['run_r'], IPHASQC['run_i'], IPHASQC['run_ha']))
fields = np.concatenate((IPHASQC['id'], IPHASQC['id'], IPHASQC['id']))
qflags = np.concatenate((IPHASQC['qflag'], IPHASQC['qflag'], IPHASQC['qflag']))
qcproblems = np.concatenate((IPHASQC['problems'], IPHASQC['problems'], IPHASQC['problems']))
depth5sig = np.concatenate((IPHASQC['r5sig_judged'],
                            IPHASQC['i5sig_judged'],
                            IPHASQC['h5sig_judged']))

field_dict = dict(zip(runs, fields))
qflag_dict = dict(zip(runs, qflags))
qcproblems_dict = dict(zip(runs, qcproblems))
depth5sig_dict = dict(zip(runs, depth5sig))

# Add the IPHAS field number
field = Column(name='fieldid', data=[field_dict[r] for r in t['run']])
t.add_column(field)

# Add the DR2 quality grade
qcgrade = Column(name='qcgrade', data=[qflag_dict[r] for r in t['run']])
t.add_column(qcgrade)

# Add the 'quality problems' summary
qcproblems = Column(name='qcproblems', data=[qcproblems_dict[r] for r in t['run']])
t.add_column(qcproblems)

# Add the 5-sigma detection limit
depth = Column(name='depth', data=[depth5sig_dict[r] for r in t['run']])
t.add_column(depth)

# Limit the number of decimals in the ascii output:
t['ra'].format = '{0:.3f}'
t['dec'].format = '{0:.3f}'

t.remove_column('airmass')
t.sort(['run', 'ccd'])

# We will export the resulting table to FITS, ASCII, and SQLITE
# First, export to FITS
columns = ['run', 'ccd', 'url', 'ra', 'dec', 'band', 'utstart',
           'fieldid', 'in_dr2', 'qcgrade', 'qcproblems',
           'exptime', 'seeing', 'elliptic', 'skylevel', 'skynoise',
           'depth', 'photzp', 'confmap',
           'ra_min', 'ra_max', 'dec_min', 'dec_max']
t[columns].write('iphas-images.fits.gz', overwrite=True)

# Export to ASCII
t['url', 'ra', 'dec', 'band', 'fieldid', 'in_dr2', 'qcgrade'].write('iphas-images.txt', format='ascii.fixed_width')

# Export to SQLITE (using atpy as astropy doesn't support sqlite yet)
import atpy
tbl = atpy.Table('iphas-images.fits.gz', name='images')
tbl.write('sqlite', 'iphas-images.sqlite', overwrite=True)

# For fast queries, you might want to do:
# CREATE INDEX images_ra_min_idx ON images(ra_min);
# CREATE INDEX images_ra_max_idx ON images(ra_max);
# VACUUM;