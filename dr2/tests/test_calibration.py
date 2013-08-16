import numpy as np
from numpy import all, abs
from .. import calibration


def test_glazebrook_equation():
    """Tests the example given in the Glazebrook paper

    (note that the example in the paper is slightly wrong)
    """
    # Setup the example from the Glazebrook paper
    runs = np.array([1, 2, 3, 4, 5, 6])
    anchors = np.array([False, False, False, False, True, True])
    overlaps = {1: {'runs': [2, 6], 'offsets': [0.5, 1.25],    'weights': [1, 1]},
                2: {'runs': [1, 6], 'offsets': [-0.5, 0.75],   'weights': [1, 1]},
                3: {'runs': [4],    'offsets': [-1.0],         'weights': [1]},
                4: {'runs': [3, 5], 'offsets': [+1.0, 1.0],    'weights': [1, 1]},
                5: {'runs': [4],    'offsets': [-1.0],         'weights': [1]},
                6: {'runs': [1, 2], 'offsets': [-1.25, -0.75], 'weights': [1, 1]}}

    class Calibration(object):
        def __init__(self, runs, overlaps, anchors):
            self.runs = runs
            self.overlaps = overlaps
            self.anchors = anchors
        def get_runs(self):
            return self.runs
        def get_overlaps(self):
            return self.overlaps
        def get_anchors(self):
            return self.anchors

    cal = Calibration(runs, overlaps, anchors)

    g = calibration.Glazebrook(cal)
    g.solve()

    A_expected = np.matrix([[-2.,  1.,  0.,  0.],
                            [ 1., -2.,  0.,  0.],
                            [ 0.,  0., -1.,  1.],
                            [ 0.,  0.,  1., -2.]])
    assert(all(g.A.todense() == A_expected))

    print(g.A.todense())

    b_expected = np.array([1.75, 0.25, -1., 2.])
    assert(all(abs(g.b - b_expected) < 1e-7))

    solution_expected = [-1.25, -0.75, 0.0, -1.0]
    assert(all(abs(g.solution[0] - solution_expected) < 1e-7))
