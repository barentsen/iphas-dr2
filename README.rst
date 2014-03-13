====================
IPHAS Data Release 2
====================

This repository contains the source code that was used to uniformly calibrate
and compile the source catalogue for the Second Data Release of the
INT/WFC Photometric H-Alpha Survey of the Northern Galactic Plane (IPHAS DR2).

The principal task of this code is to take the detection tables provided by the source
detection and measurement pipeline of the Cambridge Astronomical Survey Unit (CASU),
and transform them into a globally calibrated source catalogue which details the best
measurement of each unique source.

More information about the project can be found at www.iphas.org

Contents
--------
- ``pipeline.py``: the master script used to create the data release.
- ``paper/``: the accompanying paper describing the generated catalogue.
- ``dr2/``: the Python package that is used by `pipeline.py` to create 
  the data release, in particular the source catalogue.
  It consists of several modules:
  
  + ``dr2.detections``: module used to convert the detection tables generated 
    by the CASU pipeline into tables with user-friendly columns,
    broadly following the WSA/UKIDSS conventions.
  + ``dr2.bandmerging``: modules used to merge r/i/Halpha detection tables 
    into band-merged field catalogues.
  + ``dr2.offsets``: computes the magnitude shifts between exposure overlaps.
  + ``dr2.calibration``: computes a global photometric calibration from the
    data provided by offsets.py.
  + ``dr2.seaming``: identifies multiple detections of the same object 
    and identifies the best (primary) detection for each unique source.
  + ``dr2.concatenating``: concatenates the output of the seaming module
    into the final source catalogue.
  + ``dr2.images``: creates a release of the pipeline-processed images
    with update headers to reflect the IPHAS DR2 uniform calibration.
- ``scripts/``: various helper scripts, including:

  + ``postgresql/``: script to ingest the DR2 catalogue into a PostgreSQL db.
  + ``sqlite/``: script to ingest the DR2 catalogue into an SQLite db.

Usage
-----
The script ``pipeline.py`` in the root directory contains all instructions
necessary to produce the IPHAS DR2 source catalogue.
It depends on the IPython.Parallel package and several modules from the included ``dr2``
package to generate the data release in parallel on a computing cluster.

Dependencies
------------
IPHAS DR2 was generated using

- Python 2.7.6
- Astropy 0.3
- Numpy 1.7.1
- Scipy 0.12.0
- IPython 0.13.2


License
-------
The source code is released under the MIT license, see the ``LICENSE`` file.
