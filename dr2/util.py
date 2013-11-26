#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Utility functions for photometric data releases.

This module contains a range of utility functions needed for processing 
IPHAS/VPHAS data release products.
"""
from __future__ import division, print_function, unicode_literals
import numpy as np
import socket
import os

from dr2 import constants

__author__ = 'Geert Barentsen'
__copyright__ = 'Copyright, The Authors'
__credits__ = ['Geert Barentsen', 'Hywel Farnhill', 'Janet Drew']


def sphere_dist(lon1, lat1, lon2, lat2):
    """
    Haversine formula for angular distance on a sphere: more stable at poles.

    Inputs must be in DEGREES.
    Result is also in DEGREES.

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
    # Bugfix: crossing meridian
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

    Parameters
    ----------
    ra : float [degrees]
        Right Ascension of the position to match.

    dec : float [degrees]
        Declination of the position to match.

    ra_array : array of floats [degrees]
        Array of candidate positions.

    dec_array : array of floats [degrees]
        Array of candidate positions.

    matchdist : float [arcsec]
        Maximum matching distance.

    Returns
    -------
    idx : integer or None
        idx of the object in ra_array/dec_array which most closely matches
        the position ra/dec. If no match is found within the maximum matching 
        distance, then None is returned.
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


def get_pid():
    """Returns the hostname and process identifier.

    Returns
    -------
    pid : string
        A string of the form "hostname/process_id".
    """
    pid = '{0}/{1}'.format(socket.gethostname(),
                           os.getpid())
    return pid


def setup_dir(path):
    """Setup an output directory, i.e. make sure it exists.

    Parameters
    ----------
    path : string
        Directory to create if it does not already exist.
    """
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except OSError:  # "File already exist" can occur due to parallel running
        pass

def run2field(run, band):
    """Convert a run number to the field identifier."""
    assert(band in constants.BANDS)
    idx = np.where(constants.IPHASQC['run_{0}'.format(band)] == run)
    if len(idx) > 0:
        return constants.IPHASQC['id'][idx[0]][0]
    else:
        return None
