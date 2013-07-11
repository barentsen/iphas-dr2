#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Produces IPHAS Data Release 2 using an MPI computing cluster."""
import os
from IPython import parallel

__author__ = 'Geert Barentsen'

# Create the cluster view
client = parallel.Client(profile='mpi')
cluster = client.load_balanced_view()

# Sync imports across all nodes
with client[:].sync_imports():
    # Make sure the IPHAS DR2 module is in the path
    import sys
    sys.path.append('/home/gb/dev/iphas-dr2')
    client[:].execute("sys.path.append('/home/gb/dev/iphas-dr2')", block=True)
    # Import DR2 generation modules
    from dr2 import constants
    from dr2 import detections



detections.create_index(cluster,
                        data=os.path.join(constants.RAWDATADIR, 'iphas_sep2005'),
                        target=os.path.join(constants.DESTINATION, 'runs.csv'))
#detections.sanitise_zeropoints()         # Produces zeropoints.csv

"""
detections.create_catalogues(cluster,
                             data=os.path.join(constants.RAWDATADIR, 'iphas_sep2005'),
                             target=os.path.join(constants.DESTINATION, 'detected'))
"""

"""
detections.create_catalogues(directory)  # Single-filter catalogues
offsets.compute_offsets()
calibration.run_glazebrook()             # Re-calibration (minimises offsets)
bandmerging.bandmerge()                  # Band-merge
concatenation.concatenate()
"""
