""" 
This module implements the Alternating Regression code which is as 
described in 
"Parameter Estimation in S-systems with Alternating Regression" - 
Chou,Voit, et. al
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1586003/

Input  : 
An S-system model encoded with the type of experiment encoded
This model is mostly meant to work with the Chou2006 S-system model

"""

from utility import basedir,logdir


class ARSolver(object) :
	"""
Houses the Alternating Regression core routine and supporting scripts
to enforce constraints and good behavior of algorithm
	"""
	def __init__(self,ss=None) :  
		""" Init the AR solver""" 
		self.ss = ss
		
	def	_preprocessor(self) :  
		""" Run preprocessing on ss to make compatible with exp type
		 """
		pass

	def _enforce_cons(self) : 
		""" Enforce constraints"""
		pass
	
	def _monitor(self) : 
		""" Monitor and enforce constraints as needed """
		pass

	def solve(self) : 
		""" Runs the core routine """
		pass		

	def _core(self) : 
		""" Core routine of the ARSolver class"""
		pass


class TrajectoryTracker(object) : 
	""" Tracks trajectory of the solution path """
	def __init__(self) : 
		pass


if __name__ == '__main__' : 
	ar = ARSolver() 
