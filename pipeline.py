#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pipeline script to produce IPHAS Data Release 2 using an MPI cluster.

This script drives the entire pipeline for converting the single-band 
catalogues produced by the Cambridge Astronomical Surveys Unit (CASU) 
into a band-merged and homogeneously calibrated source catalogue 
of unique objects.

This script relies on the presence of a running computing cluster setup with 
IPython.parallel (using the "ipcluster" tool); i.e. before running the pipeline,
execute 'qsub scripts/pipeline-cluster.pbs'.
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

# While in development, reload every module by default
# This is to make sure that the latest version gets used
client[:].execute('reload(constants)', block=True)
client[:].execute('reload(util)', block=True)
client[:].execute('reload(detections)', block=True)
client[:].execute('reload(offsets)', block=True)
client[:].execute('reload(calibration)', block=True)
client[:].execute('reload(bandmerging)', block=True)
client[:].execute('reload(seaming)', block=True)
client[:].execute('reload(concatenating)', block=True)

"""
Pipeline starts here
"""
# Create an index of all single-band catalogues
# detections.create_index(cluster)  # produces 'runs.csv'

# Enforce zp(Halpha) = zp(r) - 3.14
#detections.sanitise_zeropoints()  # produces 'zeropoint-overrides.csv'

# Convert the single-band catalogues from CASU into our own catalogue format
#detections.convert_catalogues(cluster)  # produces 'detected/nnnnnnn_det.fits'

# Bandmerge all runs obtained at the same epoch and pointing
#bandmerging.bandmerge(cluster)  # produces 'bandmerged/nnnn.fits'

# Compute the magnitude offsets between all runs; necessary for re-calibration
offsets.compute_offsets(cluster)  # produces 'offsets-{r|i|ha}.csv'

# Find the set of zeropoint shifts which minimize the offsets obtained above
#calibration.calibrate()  # produces 'calibration/calibration-{r|i|ha}.csv'

# Apply the zeropoint shifts found above to the bandmerged catalogues
#calibration.apply_calibration(cluster) # produces 'bandmerged-calibrated/nnnn.fits'

# Identify duplicate detections where multiple pointings overlap ('seams');
# this requires up to 4 GB RAM per process,
# hence we only use only every fourth cluster. node (each having 1 GB).
#cluster_highmem = client[::3]
#seaming.seam(cluster_highmem)  # produces 'seamed/nnnn.fits'

# Finally, concatenate the individual pointings into a single catalogue 
#concatenating.concatenate(cluster_highmem)
