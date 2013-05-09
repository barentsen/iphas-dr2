"""
Input: list of runs with overlaps and initial zeropoints
Output: zeropoints

Input
=====

"""
import numpy as np
from scipy import sparse
from scipy.sparse import linalg
from astropy.io import fits, ascii
from astropy import log


class Glazebrook(object):
    """Fit a global calibration using the method by (Glazebrook et al. 1994)

    Uses sparse matrix functions (scipy.sparse) for efficiency.
    """

    def __init__(self, runs, overlaps, anchors):
        self.runs = runs
        self.overlaps = overlaps
        self.anchors = anchors

        self.n_nonanchors = (~anchors).sum()

    def _A(self):
        log.info('Creating a sparse {0}x{0} matrix'.format(self.n_nonanchors))
        A = sparse.lil_matrix((self.n_nonanchors,
                               self.n_nonanchors))
        for i, run in enumerate(self.runs[~self.anchors]):
            for j, run2 in enumerate(self.runs[~self.anchors]):
                if j < i:  # The matrix is symmetric
                    continue
                if i == j:
                    A[i, j] = -len(self.overlaps[run2]['runs'])
                elif run2 in self.overlaps[run]['runs']:
                    A[i, j] = 1
                    A[j, i] = 1
        return A

    def _b(self):
        b = np.zeros(self.n_nonanchors)
        for i, run in enumerate(self.runs[~self.anchors]):
            b[i] = np.sum(self.overlaps[run]['offsets'])
        return b

    def solve(self):
        self.A = self._A()
        self.b = self._b()
        log.info('Now solving the matrix equation')
        # Note: there should be alternative algorithms for symmetric
        # matrices which are faster.
        self.solution = linalg.lsqr(self.A, self.b,
                                    atol=1e-12, iter_lim=3e5)
        return self.solution

    def write(self, filename):
        f = open(filename, 'w')
        f.write('run,shift\n')
        for i, myrun in enumerate(self.runs[~self.anchors]):
            f.write('{0},{1}\n'.format(myrun, self.solution[0][i]))
        for myrun in self.runs[anchors]:
            f.write('{0},{1}\n'.format(myrun, 0.0))
        f.close()


def prepare_glazebrook_data(band='r'):
    assert(band in ['r', 'i', 'ha'])
    IPHASQC = fits.getdata('/home/gb/dev/iphas-qc/qcdata/iphas-qc.fits', 1)

    # Parse data
    filename_offsets = 'glazebrook/offsets-{0}.csv'.format(band)
    log.info('Reading {0}'.format(filename_offsets))
    offsetdata = ascii.read(filename_offsets)

    # All the runs, sorted numerically
    runs = np.sort(np.unique(offsetdata['run1']))

    # 'anchors' is a boolean array indicating anchor status
    anchors = []
    QC_ANCHORS = IPHASQC.field('anchor')
    QC_RUNS = IPHASQC.field('run_{0}'.format(band))
    for myrun in runs:
        if QC_ANCHORS[QC_RUNS == myrun][0] == 1:
            anchors.append(True)
        else:
            anchors.append(False)
    anchors = np.array(anchors)

    # Dictionary of field overlaps
    overlaps = {}
    for row in offsetdata:
        myrun = row['run1']
        if myrun not in overlaps:
            overlaps[myrun] = {'runs': [], 'offsets': []}
        overlaps[myrun]['runs'].append( row['run2'] )
        overlaps[myrun]['offsets'].append( row['offset'] )

    return (runs, overlaps, anchors)


if __name__ == '__main__':
    band = 'r'
    runs, overlaps, anchors = prepare_glazebrook_data(band)
    g = Glazebrook(runs, overlaps, anchors)
    s = g.solve()[0]
    g.write('glazebrook/calibration-{0}.csv'.format(band))
