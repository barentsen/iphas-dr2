IPHAS Data Release 2
====================

Python package, scripts and documentation used to produce IPHAS Data Release 2.

Contents
--------
* dr2: Python DR2 package, consisting of the following modules:
  * detections.py: converts CASU pipeline output to tables with IPHAS columns
  * offsets.py: computes the photometric shifts between exposure overlaps
  * calibration.py: computes a global photometric calibration from the offsets
  * bandmerging.py: merges the H-alpha/r/i detections for the same pointing
  * seaming.py: identifies multiple detections of the same source
  * concatenating.py: produces the final Primary Source Catalogue
* scripts: scripts used to run the pipeline on a computing cluster.
* documentation: description of tables.

Workflow
--------
To carry out a global photometric calibration, run the following modules:
detections.py => offsets.py => calibration.py

To generate the catalogues:
detections.py => bandmerging.py => seaming.py => concatenating.py

Dependencies
------------
astropy, numpy, scipy.

License
--------
Copyright, the authors 2013.
You are required to contact the authors if you wish to use the contents of this repository.
