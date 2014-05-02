#!/usr/bin/env python
"""Copies the confidence maps into the data release directory."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
__author__ = "Geert Barentsen"

import os
import subprocess
import numpy as np
from astropy import table
from astropy import log
from dr2 import constants

if __name__ == '__main__':
    destination = os.path.join(constants.PATH_IMAGES, 'confmaps')
    imagetable = table.Table.read(os.path.join(constants.PATH_IMAGES,
                                           'iphas-images.fits'))

    for filename in np.unique(imagetable['confmap']):
        source = os.path.join(constants.RAWDATADIR, filename)
        target = os.path.join(destination, filename)
        if not os.path.exists(os.path.dirname(target)):
            log.info('Creating {0}'.format(os.path.dirname(target)))
            os.makedirs(os.path.dirname(target))
        log.info('Copying {0} to {1}'.format(source, target))
        c = subprocess.call(['cp', '-a', source, target])
