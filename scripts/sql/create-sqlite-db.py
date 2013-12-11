"""
Create and query SQLite files from survey data

"""
from astropy.io import fits
import sqlite3
import os
import glob
import numpy as np
from astropy import log

DATADIR = '/car-data/gb/iphas-dr2-rc4/concatenated/full'
DESTINATION = '/car-data/gb/iphas-dr2-rc4'

class SurveyDB(object):

    def __init__(self, filename):
        # sqlite doesn't know how to handle certain numpy types unless told
        sqlite3.register_adapter(np.int8, int)
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

    def query(self, sql):
        #self.cursor.execute(sql)
        #astropy.table.Table(np.array(self.cursor.fetchall())).write('test.fits')
        pass

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
        log.info('Indexing (ra,dec)')
        self.cursor.execute("PRAGMA temp_store=FILE")
        self.cursor.execute("PRAGMA temp_store_directoryi='/car-data/gb/tmp'")
        self.cursor.execute("PRAGMA cache_size = '2000000'")
        self.cursor.execute('CREATE INDEX iphas_ra_dec_idx ON iphas(ra, dec)')
        self.cursor.execute('CREATE INDEX iphas_l_b_idx ON iphas(l, b)')

    def insert_ndarray(self, data, columns):
        log.info('Ingesting data')
        self.cursor.executemany('INSERT INTO iphas VALUES (' 
                                + ','.join(['?']*len(columns)) + ')',
                                data)        

    def insert_fits(self, filename, columns):
        log.info('Ingesting {0}'.format(filename))
        # Using .__array__() is important for speed
        self.cursor.executemany('INSERT INTO iphas VALUES (' 
                                + ','.join(['?']*len(columns)) + ')',
                                fits.getdata(filename, 1).__array__())



def create_iphas_light():
    """Running cost: 91m22.244s (6 Nov 2013)"""
    cols = ['name', 'ra', 'dec', 
               'r', 'rErr', 'i', 'iErr', 'ha', 'haErr',
               'mergedClass', 'errBits']
    db = SurveyDB('iphas-dr2-light.db')
    db.optimise_inserts()
    db.create_table('iphas', cols)
    for filename in np.sort(glob.glob('/home/gb/tmp/iphas-dr2-rc3/concatenated/light/*-light.fits.gz')):
        db.insert_fits(filename, cols)
    db.create_indexes()
    db.commit()

def create_iphas_full():
    """Running cost: about 10 hours?"""
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
    db = SurveyDB(os.path.join(DESTINATION, 'iphas-dr2-full.db'))
    db.optimise_inserts()
    """
    db.create_table('iphas', cols)
    for filename in np.sort(glob.glob(os.path.join(DATADIR, '*.fits.gz'))):
        db.insert_fits(filename, cols)
    db.commit()
    """
    db.create_indexes()
    db.commit()
    

if __name__ == '__main__':
    #create_iphas_light()
    create_iphas_full()

#data = np.asarray(fits.getdata(mycatalogue, 1))
#data = fits.getdata(mycatalogue, 1).__array__()
