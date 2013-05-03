"""
Input: list of runs with overlaps and initial zeropoints
Output: zeropoints

Input
=====

"""
import numpy as np
from scipy import sparse
from scipy.sparse import linalg



class Glazebrook(object):

	def __init__(self, runs, overlaps, anchors):
		self.runs = runs
		self.overlaps = overlaps
		self.anchors = anchors

	def _A(self):
		A = sparse.lil_matrix((len(self.runs), 
			                        len(self.runs)))

		for i, run in enumerate(self.runs):
			for j, run2 in enumerate(self.runs):
				if self.anchors[i] and i == j:
					A[i,j] = 1
				elif self.anchors[i]:
					continue
				elif i == j:
					A[i,j] = 1 - (len(self.overlaps[run2]['runs']) + 1)
				elif run2 in self.overlaps[run]['runs']:
					A[i,j] = 1

		return A

	def _b(self):
		b = np.zeros(len(self.runs))
		for i, run in enumerate(self.runs):
			if not self.anchors[i]:
				b[i] = np.sum(self.overlaps[run]['offsets'])
		return b

	def solve(self):
		A = self._A()
		b = self._b()
		s = linalg.lsqr(A, b)
        return s




if __name__ == '__main__':
	runs = np.array([1, 2, 3, 4, 5, 6])
	anchors = np.array([False, False, False, False, True, True])
	overlaps = {1: {'runs':[2, 6], 'offsets':[0.5, 1.25]},
	            2: {'runs':[1, 6], 'offsets':[-0.5, 0.75]},
	            3: {'runs':[4], 'offsets':[-1.0]},
	            4: {'runs':[3, 5], 'offsets':[+1.0, 1.0]},
	            5: {'runs':[4], 'offsets':[-1.0]},
	            6: {'runs':[1,2], 'offsets':[-1.25,-0.75]}}

	g = Glazebrook(runs, overlaps, anchors)
	A = g._A()
	b = g._b()
	s = g.solve()[0]