#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Produces IPHAS Data Release 2 using an MPI computing cluster."""
from IPython import parallel

__author__ = 'Geert Barentsen'

# Create the cluster view
client = parallel.Client('/home/gb/.config/ipython/profile_mpi/security/ipcontroller-pipeline-client.json')
cluster = client.load_balanced_view()

# Sync imports across all nodes
with client[:].sync_imports():
    # Make sure the IPHAS DR2 module is in the path
    import os
    import sys
    sys.path.append('/home/gb/dev/iphas-dr2')
    client[:].execute("sys.path.append('/home/gb/dev/iphas-dr2')", block=True)
    # Import DR2 generation modules
    from dr2 import constants
    from dr2 import detections
    from dr2 import offsets
    from dr2 import bandmerging



#detections.create_index(cluster)
#data=os.path.join(constants.RAWDATADIR, 'iphas_sep2005'),
#detections.sanitise_zeropoints()         # Produces zeropoint-overrides.csv
detections.convert_catalogues(cluster)
#offsets.compute_offsets(cluster)
#bandmerging.bandmerge(cluster)



"""
calibration.run_glazebrook()             # Re-calibration (minimises offsets)
concatenation.concatenate()
"""
