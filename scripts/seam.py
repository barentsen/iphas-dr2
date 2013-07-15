#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Produces IPHAS Data Release 2 using an MPI computing cluster."""
from IPython import parallel
from astropy import log

__author__ = 'Geert Barentsen'

# Create the cluster view
client = parallel.Client('/home/gb/.config/ipython/profile_mpi/security/ipcontroller-pipeline-client.json')
cluster = client[:]
log.info('Using {0} cores'.format(len(cluster)))

# Sync imports across all nodes
with client[:].sync_imports():
    # Make sure the IPHAS DR2 module is in the path
    import os
    import sys
    sys.path.append('/home/gb/dev/iphas-dr2')
    client[:].execute("sys.path.append('/home/gb/dev/iphas-dr2')", block=True)
    # Import DR2 generation modules
    from dr2 import constants
    from dr2 import seaming

client[:].execute('reload(constants)', block=True)
client[:].execute('reload(seaming)', block=True)

seaming.seam(cluster)
