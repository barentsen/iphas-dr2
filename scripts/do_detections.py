import sys
import os
from IPython import parallel
from astropy import log

sys.path.append('/home/gb/dev/iphas-dr2')
from dr2 import constants
from dr2 import detections

log.setLevel('DEBUG')

client = parallel.Client(profile='mpi')
cluster = client.load_balanced_view()


def create_index(clusterview,
                 target=os.path.join(constants.DESTINATION, 'runs.csv'),
                 data=constants.RAWDATADIR):
    """Produces a CSV file detailing the properties of all runs."""
    out = detections.index_setup(target)
    catalogues = detections.list_catalogues(data)

    # Index each pipeline catalogue
    results = clusterview.imap(detections.index_one, catalogues)  # returns an iterator
    for r in results:
        if r is None:
            continue
        out.write(r+'\n')
        out.flush()
    out.close()


create_index(cluster,
             data=os.path.join(constants.RAWDATADIR, 'iphas_sep2005'),
             target=os.path.join(constants.DESTINATION, 'runs.csv'))
#detections.sanitise_zeropoints()         # Produces zeropoints.csv
#detections.create_catalogues(cluster,
#                             target=os.path.join(constants.DESTINATION, 'detected'))

"""
detections.create_catalogues(directory)  # Single-filter catalogues
offsets.compute_offsets()
calibration.run_glazebrook()             # Re-calibration (minimises offsets)
bandmerging.bandmerge()                  # Band-merge
concatenation.concatenate()
"""
