import sqlite3
from astropy import log
import numpy as np
import math
import os
from .. import util
from .. import constants


# Numpy dtype for a row in the IPHAS DR2 full source catalogue
mydtype = [('name', 'S19'), ('ra', '>f8'), ('dec', '>f8'), ('sourceID', 'S14'), 
           ('posErr', '>f4'), ('l', '>f8'), ('b', '>f8'), 
           ('mergedClass', '>i2'), ('mergedClassStat', '>f4'), ('pStar', '>f4'), 
           ('pGalaxy', '>f4'), ('pNoise', '>f4'), ('pSaturated', '>f4'), 
           ('rmi', '>f4'), ('rmha', '>f4'), 
           ('r', '>f4'), ('rErr', '>f4'), 
           ('rPeakMag', '>f4'), ('rPeakMagErr', '>f4'), 
           ('rAperMag1', '>f4'), ('rAperMag1Err', '>f4'), 
           ('rAperMag3', '>f4'), ('rAperMag3Err', '>f4'), 
           ('rGauSig', '>f4'), ('rEll', '>f4'), ('rPA', '>f4'), 
           ('rClass', '>i2'), ('rClassStat', '>f4'), 
           ('rErrBits', '>i2'), ('rMJD', '>f8'), ('rSeeing', '>f4'), 
           ('rDetectionID', 'S15'), ('rX', '>f4'), ('rY', '>f4'), 
           ('i', '>f4'), ('iErr', '>f4'), 
           ('iPeakMag', '>f4'), ('iPeakMagErr', '>f4'), 
           ('iAperMag1', '>f4'), ('iAperMag1Err', '>f4'), 
           ('iAperMag3', '>f4'), ('iAperMag3Err', '>f4'), 
           ('iGauSig', '>f4'), ('iEll', '>f4'), ('iPA', '>f4'), 
           ('iClass', '>i2'), ('iClassStat', '>f4'), 
           ('iErrBits', '>i2'), ('iMJD', '>f8'), ('iSeeing', '>f4'), 
           ('iDetectionID', 'S15'), ('iX', '>f4'), ('iY', '>f4'), 
           ('iXi', '>f4'), ('iEta', '>f4'), 
           ('ha', '>f4'), ('haErr', '>f4'), 
           ('haPeakMag', '>f4'), ('haPeakMagErr', '>f4'), 
           ('haAperMag1', '>f4'), ('haAperMag1Err', '>f4'), 
           ('haAperMag3', '>f4'), ('haAperMag3Err', '>f4'), 
           ('haGauSig', '>f4'), ('haEll', '>f4'), ('haPA', '>f4'), 
           ('haClass', '>i2'), ('haClassStat', '>f4'), 
           ('haErrBits', '>i2'), ('haMJD', '>f8'), ('haSeeing', '>f4'), 
           ('haDetectionID', 'S15'), ('haX', '>f4'), ('haY', '>f4'), 
           ('haXi', '>f4'), ('haEta', '>f4'), 
           ('brightNeighb', 'i1'), ('deblend', 'i1'), ('saturated', 'i1'), 
           ('errBits', '>i2'), ('nBands', '>i2'), ('reliable', 'i1'), 
           ('fieldID', 'S14'), ('fieldGrade', 'S3'), ('night', '>i4'), 
           ('seeing', '>f4'), ('ccd', '>i2'), ('nObs', '>i2'), 
           ('sourceID2', 'S15'), ('fieldID2', 'S15'), 
           ('r2', '>f4'), ('rErr2', '>f4'), ('i2', '>f4'), ('iErr2', '>f4'), 
           ('ha2', '>f4'), ('haErr2', '>f4'), ('errBits2', '>i4')]



class SourceCatalogueDB(object):

    def __init__(self):
        #self.filename = os.path.join(constants.PACKAGEDIR, 'io', 'iphas-dr2-full.db')
        self.filename = os.path.join(constants.DESTINATION, 'iphas-dr2-full.db')
        self.connection = sqlite3.connect(self.filename)
        self.cursor = self.connection.cursor()

    def query(self, sql):
        """Returns numpy ndarray."""
        return np.fromiter(self.cursor.execute(sql), mydtype)

    def conesearch(self, ra, dec, radius=60.0):
        """
        Radius in arcsec
        """
        sr_degrees = radius/3600.
        sql = """
                SELECT * FROM iphas 
                WHERE ra BETWEEN %s-%s AND %s+%s 
                AND dec BETWEEN %s-%s AND %s+%s;
        """ % (ra, sr_degrees*1.1*math.cos(math.radians(dec)), ra, sr_degrees*1.1*math.cos(math.radians(dec)),
               dec, sr_degrees, dec, sr_degrees,)
        log.info(sql)
        result = self.query(sql)
        sep = util.sphere_dist(ra, dec, result['ra'], result['dec'])
        return result[sep < sr_degrees]

    def boxsearch(self, lon1, lon2, lat1, lat2, cond=''):
        """
        """
        sql = """
                SELECT * FROM iphas 
                WHERE l >= %s AND l <= %s 
                AND b >= %s AND b <= %s
                %s;""" % (lon1, lon2, lat1, lat2, cond)
        log.info(sql)
        return self.query(sql)


if __name__ == '__main__':
    #r = conesearch(ra, dec)
    #r = box(29.1, 30.1, -0.1, 0.1)
    r = plotboxes()