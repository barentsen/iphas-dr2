"""
Input: list of runs with overlaps and initial zeropoints
Output: zeropoints

Input
=====

"""
import numpy as np
from scipy import sparse



class Glazebrook(object):

	def __init__(self, runs):
		self.runs = runs

	def _create_A(self):
		self.A = sparse.lil_matrix((len(self.runs), 
			                        len(self.runs)))



if __name__ == '__main__':
	runs = [{'run':1, 'zp':20.0},
	        {'run':2, 'zp':20.0}]

	g = Glazebrook(runs)
	g._create_A()