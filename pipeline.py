#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pipeline script used to produce IPHAS Data Release 2 on a computing cluster.

Summary
-------
This script drives the entire creation of the IPHAS DR2 source catalogue.
Its purpose is to convert the IPHAS single-band detection tables, which are
kindly provided by the Cambridge Astronomical Surveys Unit (CASU),
into a band-merged and globally calibrated source catalogue which lists only
the best-available measurement for each unique source.

Usage
-----
This script relies on the presence of a running instance of an IPython.parallel
cluster, which can be started on a local machine using:
    
    ``ipcluster start``

or can be started on a PBS-controlled computing cluster by submitting the job:

   ``qsub scripts/pipeline-cluster.pbs``.

Once an IPython.parallel cluster is running, you can simply run

   ``python pipeline.py``

...and the entire data release will be created within a few thousand CPU hours.
"""
from IPython import parallel
from astropy import log

__author__ = 'Geert Barentsen'

# Create the cluster view
MYCLUSTER = '/home/gb/.config/ipython/profile_mpi/security/ipcontroller-pipeline-client.json'
client = parallel.Client(MYCLUSTER)

# We use a direct view because the load_balanced_view() was found to be buggy
cluster = client[:]
log.info('Using {0} cores'.format(len(cluster)))

# Some tasks require up to 4 GB RAM per process,
# we will execute them on every fourth core only 
# as each core is requested to have 1 GB available.
# This is a hack to be able to do everything on one cluster job reservation.
cluster_highmem = client[::4]
log.info('Using {0} high-memory cores'.format(len(cluster_highmem)))

# Import our required modules across all computing nodes
with client[:].sync_imports():
    import os
    import sys
    # Make sure the IPHAS DR2 module is in the path at all nodes
    sys.path.append('/home/gb/dev/iphas-dr2')
    client[:].execute("sys.path.append('/home/gb/dev/iphas-dr2')", block=True)
    from dr2 import constants
    from dr2 import util
    from dr2 import detections
    from dr2 import offsets
    from dr2 import calibration
    from dr2 import bandmerging
    from dr2 import seaming
    from dr2 import concatenating
    from dr2 import images

# While in development, reload every module by default,
# to make sure that the latest version gets used. If the pipeline were ever
# to be used in production (never), this could be removed. 
client[:].execute('reload(constants)', block=True)
client[:].execute('reload(util)', block=True)
client[:].execute('reload(detections)', block=True)
client[:].execute('reload(offsets)', block=True)
client[:].execute('reload(calibration)', block=True)
client[:].execute('reload(bandmerging)', block=True)
client[:].execute('reload(seaming)', block=True)
client[:].execute('reload(concatenating)', block=True)
client[:].execute('reload(images)', block=True)

"""
Pipeline starts here
"""
# Create an index of all single-band catalogues
detections.save_metadata(cluster)  # produces 'metadata.fits'

# Enforce zp(Halpha) = zp(r) - 3.14
detections.sanitise_zeropoints()  # produces 'zeropoints-precalibration.csv'

# Convert the single-band catalogues from CASU into our own catalogue format
detections.convert_catalogues(cluster)  # produces 'detected/nnnnnnn_det.fits'

# Bandmerge all runs obtained at the same epoch and pointing
bandmerging.bandmerge(cluster)  # produces 'bandmerged/nnnn.fits'

# Compute the magnitude offsets between all runs, which is a necessary
# input to the re-calibration step. Executing this on too many cores has been 
# found to result in # "[Errno 105] No buffer space available", 
# so we run it on cluster_highmem defined earlier.
offsets.compute_offsets(cluster_highmem)  # produces 'offsets-{r|i|ha}.csv'

# Find the set of zeropoint shifts which minimize the offsets obtained above
calibration.calibrate()  # produces 'calibration/calibration-{r|i|ha}.csv'

# Apply the zeropoint shifts found above to the bandmerged catalogues
calibration.apply_calibration(cluster) # produces 'bandmerged-calibrated/nnnn.fits'

# Identify duplicate detections where multiple pointings overlap ('seams');
seaming.seam(cluster_highmem)  # produces 'seamed/nnnn.fits'

# Finally, concatenate the individual pointings into a single catalogue 
concatenating.concatenate(cluster_highmem)

# Prepare images for release
images.prepare_images(cluster)
