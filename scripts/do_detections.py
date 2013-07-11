sys.path.append('/home/gb/dev/iphas-dr2')
from dr2 import constants
from dr2 import detections
from IPython import parallel

client = parallel.Client(profile='mpi')
clusterview = client.load_balanced_view()

detections.index_all_parallel(constants.RAWDATADIR, clusterview)
#detections.convert_all_parallel(clusterview, constants.RAWDATADIR)

"""
detections.sanitise_zeropoints()         # Produces zeropoints.csv
detections.create_catalogues(directory)  # Single-filter catalogues
offsets.compute_offsets()
calibration.run_glazebrook()             # Re-calibration (minimises offsets)
bandmerging.bandmerge()                  # Band-merge
concatenation.concatenate()
"""
