====================
IPHAS Data Release 2
====================

Python package, scripts and documentation used to produce IPHAS Data Release 2.

:authors: Geert Barentsen; Hywel Farnhill
:dependencies: astropy, numpy, scipy, IPython.

Contents
--------
- **dr2**: Python package consisting of the following modules:

  + **detections.py**: converts the detection catalogues from the CASU pipeline into user-friendly tables.
  + **offsets.py**: computes the photometric shifts between exposure overlaps.
  + **calibration.py**: computes a global photometric calibration from the offsets.
  + **bandmerging.py**: merges the H-alpha/r/i detections for the same pointing.
  + **seaming.py**: identifies multiple detections of the same source, and assigns the best ('primary') detection.
  + **concatenating.py**: produces the final Primary Source Catalogue.
- **scripts**: various helper scripts.
- **website**: redesigned website to be launched with DR2.

Workflow
--------
The script **pipeline.py** in the root directory contains all the commands necessary to produce the IPHAS DR2 source catalogue.

License
--------
Copyright, the authors 2013.
You are required to contact the authors if you wish to use the contents of this repository.