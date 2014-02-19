"""
Creates an SQLite database from the collection of IPHAS FITS catalogue.

"""
from astropy.io import fits
import sqlite3
import os
import glob
import numpy as np
from astropy import log
from dr2 import constants


class SurveyDB(object):

    def __init__(self, filename):
        # sqlite doesn't know how to handle certain numpy types unless told
        sqlite3.register_adapter(np.int8, bool)
        sqlite3.register_adapter(np.int16, int)
        sqlite3.register_adapter(np.int32, int)
        sqlite3.register_adapter(np.float32, float)
        #sqlite3.register_adapter(bool, int)
        #sqlite3.register_converter("BOOLEAN", lambda v: bool(int(v)))
        
        self.filename = filename
        self.connection = sqlite3.connect(self.filename,
                                          isolation_level='EXCLUSIVE')
        self.cursor = self.connection.cursor()

    def __del__(self):
        self.connection.close()

    def commit(self):
        self.connection.commit()

    def optimise_inserts(self):
        """Optimise the cursor for bulk inserts
        Inspired by http://blog.quibb.org/2010/08/fast-bulk-inserts-into-sqlite/
        """
        self.cursor.execute('PRAGMA synchronous=OFF')
        self.cursor.execute('PRAGMA count_changes=OFF')
        self.cursor.execute('PRAGMA journal_mode=MEMORY')
        self.cursor.execute('PRAGMA temp_store=MEMORY')
        self.cursor.execute('PRAGMA locking_mode=EXCLUSIVE')

    def create_table(self, name, columns):
        coldef = "('"+"', '".join(columns)+"')"
        self.cursor.execute("CREATE TABLE IF NOT EXISTS {0} {1}".format(name, coldef))

    def create_indexes(self):
        self.cursor.execute("PRAGMA temp_store=FILE")
        self.cursor.execute("PRAGMA temp_store_directory='{0}'".format(constants.TMPDIR))
        self.cursor.execute("PRAGMA cache_size = '2000000'")
        log.info('Now indexing (ra, dec)')
        self.cursor.execute('CREATE INDEX iphas_ra_dec_idx ON iphas(ra, dec)')
        log.info('Now indexing (l, b)')
        self.cursor.execute('CREATE INDEX iphas_l_b_idx ON iphas(l, b)')

    def insert_ndarray(self, data, columns):
        log.info('Ingesting data')
        self.cursor.executemany('INSERT INTO iphas VALUES (' 
                                + ','.join(['?']*len(columns)) + ')',
                                data)        

    def insert_fits(self, filename, columns):
        log.info('Ingesting {0}'.format(filename))
        # Using .__array__() is important for speed
        data = fits.getdata(filename, 1).__array__()
        # Correct boolean columns for the silly way in which FITS stores bools
        for field in ['brightNeighb', 'deblend', 'saturated', 'reliable']:
            data[field][data[field] == ord('T')] = 1
            data[field][data[field] == ord('F')] = 0
        # Insert
        self.cursor.executemany('INSERT INTO iphas VALUES (' 
                                + ','.join(['?']*len(columns)) + ')',
                                data)



def create_iphas_light():
    """Running cost: 91m22.244s (6 Nov 2013)"""
    cols = ['name', 'ra', 'dec', 
               'r', 'rErr', 'i', 'iErr', 'ha', 'haErr',
               'mergedClass', 'errBits']
    db = SurveyDB('iphas-dr2-light.db')
    db.optimise_inserts()
    db.create_table('iphas', cols)
    for filename in np.sort(glob.glob(os.path.join(constants.PATH_CONCATENATED, 'light', '*-light.fits'))):
        db.insert_fits(filename, cols)
    db.create_indexes()
    db.commit()

def create_iphas_full():
    """Creates iphas-dr2-full.db

    Running cost: about 10 hours?
    """
    cols = ['name', 'ra', 'dec',
            'sourceID', 'posErr', 'l', 'b',
            'mergedClass', 'mergedClassStat', 
            'pStar', 'pGalaxy', 'pNoise',
            'pSaturated', 'rmi', 'rmha',
            'r', 'rErr', 'rPeakMag', 'rPeakMagErr',
            'rAperMag1', 'rAperMag1Err', 'rAperMag3', 'rAperMag3Err',
            'rGauSig', 'rEll', 'rPA', 'rClass', 'rClassStat',
            'rErrBits', 'rMJD', 'rSeeing', 'rDetectionID',
            'rX', 'rY',
            'i', 'iErr', 'iPeakMag', 'iPeakMagErr',
            'iAperMag1', 'iAperMag1Err', 'iAperMag3', 'iAperMag3Err',
            'iGauSig', 'iEll', 'iPA', 'iClass', 'iClassStat',
            'iErrBits', 'iMJD', 'iSeeing', 'iDetectionID',
            'iX', 'iY', 'iXi', 'iEta',
            'ha', 'haErr', 'haPeakMag', 'haPeakMagErr',
            'haAperMag1', 'haAperMag1Err', 'haAperMag3', 'haAperMag3Err',
            'haGauSig', 'haEll', 'haPA', 'haClass', 'haClassStat',
            'haErrBits', 'haMJD', 'haSeeing', 'haDetectionID',
            'haX', 'haY', 'haXi', 'haEta',
            'brightNeighb', 'deblend', 'saturated', 'errBits',
            'nBands', 'reliable', 'fieldID', 'fieldGrade',
            'night', 'seeing', 'ccd', 'nObs', 'sourceID2',
            'fieldID2', 'r2', 'rErr2', 'i2', 'iErr2', 'ha2', 'haErr2',
            'errBits2']
    #db = SurveyDB('/home/gb/tmp/test.db')
    db = SurveyDB(os.path.join(constants.DESTINATION, 'iphas-dr2-full.db'))
    db.optimise_inserts()
    db.create_table('iphas', cols)
    files_to_insert = glob.glob(os.path.join(constants.PATH_CONCATENATED, 'full', '*.fits'))
    #files_to_insert = glob.glob(os.path.join(constants.PATH_CONCATENATED, 'full', '*215b.fits'))
    for filename in np.sort(files_to_insert):
        db.insert_fits(filename, cols)
    db.commit()
    db.create_indexes()
    db.commit()


if __name__ == '__main__':
    #create_iphas_light()
    create_iphas_full()

