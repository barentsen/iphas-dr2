#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Utility functions for data releases."""
from __future__ import division, print_function, unicode_literals
import numpy as np


def sphere_dist(lon1, lat1, lon2, lat2):
    """
    Haversine formula for angular distance on a sphere: more stable at poles.

    Inputs must be in DEGREES.

    Credit: https://github.com/astropy/astropy/pull/881
    (pull request wasn't merged at the time of writing this code,
    hence the function was copied here)
    """
    sdlat = np.sin(np.radians(lat2 - lat1) * 0.5)
    sdlon = np.sin(np.radians(lon2 - lon1) * 0.5)
    coslats = np.cos(np.radians(lat1)) * np.cos(np.radians(lat2))

    return np.degrees(2 * np.arcsin((sdlat**2 + coslats * sdlon**2) ** 0.5))


def sphere_dist_fast(lon1, lat1, lon2, lat2):
    """
    Euclidean angular distance "on a sphere" - only valid on sphere in the
    small-angle approximation.
    """

    if isinstance(lon1, np.ndarray) and len(lon1) > 1:
        lon1[(lon1 - lon2) > 180] -= 360
    elif isinstance(lon2, np.ndarray) and len(lon2) > 1:
        lon2[lon1 - lon2 > 180] -= 360
    elif lon1 - lon2 > 180:
        lon1 -= 360

    dlat = lat2 - lat1
    dlon = (lon2 - lon1) * np.cos(np.radians(lat1 + lat2) * 0.5)

    return (dlat ** 2 + dlon ** 2) ** 0.5


def crossmatch(ra, dec, ra_array, dec_array, matchdist=0.5):
    """Returns the index of the matched source.

    ra/dec in degrees.
    matching distance in arcseconds.
    """
    precut = np.abs(dec - dec_array) < (matchdist/3600.)  # Optimized
    if not precut.any():
        return None
    dist = sphere_dist(ra, dec, ra_array[precut], dec_array[precut])
    idx_closest = dist.argmin()
    if dist[idx_closest] < (matchdist / 3600.):
        return np.argwhere(precut)[idx_closest]
    else:
        return None
