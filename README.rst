====================
IPHAS Data Release 2
====================

This repository contains the source code that was used to produce the
source catalogue and the accompanying paper for the Second Data Release 
of the INT/WFC Photometric H-Alpha Survey of the Northern Galactic Plane
(IPHAS DR2).

Contents
--------
- pipeline.py: the script that steers the entire source catalogue generation
pipeline, to be used on a HPC cluster with IPython.Parallel.
- dr2/: the Python package that is used by `pipeline.py` to create 
the source catalogue. It consists of several modules:
  + detections.py: converts the detection tables generated 
  by the CASU pipeline into tables with user-friendly columns
  (broadly following UKIDSS conventions).
  + bandmerging.py: merges the r/i/Halpha detection tables into band-merged
  field catalogues.
  + offsets.py: computes the magnitude shifts between exposure overlaps.
  + calibration.py: computes a global photometric calibration from the
  data provided by offsets.py.
  + seaming.py: identifies multiple detections of the same object 
  and identifies the best ('primary') detection for each unique source.
  + concatenating.py: concatenates the output of the seaming module
  into the final source catalogue.
  + images.py: creates a release of the pipeline-processed images
  with update headers to reflect the IPHAS DR2 uniform calibration.
- scripts/: various helper scripts.
- paper/: the accompanying paper describing the data release.

Workflow
--------
The script `pipeline.py` in the root directory contains all instructions
that are necessary to produce the IPHAS DR2 source catalogue.
It uses the IPython.Parallel package and several modules from our ``dr2``
package to generate the data release in parallel on an HPC cluster.
