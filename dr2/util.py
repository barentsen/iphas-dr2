#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Utility functions."""
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
