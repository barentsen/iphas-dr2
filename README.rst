====================
IPHAS Data Release 2
====================

:summary: Python package, scripts and documentation used to produce IPHAS Data Release 2.
:authors: Geert Barentsen; Hywel Farnhill
:dependencies: astropy, numpy, scipy.

Contents
--------
- **dr2**: Python package consisting of the following modules:

  + **detections.py**: converts the detection catalogues from the CASU pipeline into user-friendly tables.
  + **offsets.py**: computes the photometric shifts between exposure overlaps.
  + **calibration.py**: computes a global photometric calibration from the offsets.
  + **bandmerging.py**: merges the H-alpha/r/i detections for the same pointing.
  + **seaming.py**: identifies multiple detections of the same source, and assigns the best ('primary') detection.
  + **concatenating.py**: produces the final Primary Source Catalogue.
- **scripts**: scripts used to run the pipeline on a computing cluster.
- **documentation**

Workflow
--------
- To carry out a global photometric calibration, call the __main__ functions of the modules in this order:
   detections.py => offsets.py => calibration.py
- To generate the IPHAS Primary Source Catalogue:
   detections.py => bandmerging.py => seaming.py => concatenating.py

License
--------
Copyright, the authors 2013.
You are required to contact the authors if you wish to use the contents of this repository.