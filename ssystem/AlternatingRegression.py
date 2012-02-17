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

from parsermanager import ParserManager
from utility import basedir,logdir
import logging, copy
from modifiers import ModifierChou2006

class ARSolver(object) :
	"""
Houses the Alternating Regression core routine and supporting scripts
to enforce constraints and good behavior of algorithm
	"""
	def __init__(self,ss=None) :  
		""" Init the AR solver""" 
		self.ss = ss
		self.logger = logging.getLogger('ss.ar')

	def _parse_initbound(self) : 
		""" Parse the soft constraints """
		logging.debug("Parsing initbound soft constraints")
		

	def _parse_initsol(self) : 
		""" Parse the initial solution """
		logging.debug("Parsing initsol initial solution")

	def _parse_modelspace(self) :
		""" Parse the modelspace which are the hard constraints """
		logging.debug("Parsing modelspace hard constraints")	

	def	_preprocessor(self) :  
		""" Run preprocessing on ss to make compatible with exp type
		 """
		logging.debug('Beginning preprocessor')
		
		# Parse entries from ss class
		self._parse_initsol()
		self._parse_initbound()
		self._parse_modelspace()
		

	def _enforce_cons(self) : 
		""" Enforce constraints"""
		logging.debug('Enforcing constraints')
		
	
	def _monitor(self) : 
		""" Monitor and enforce constraints as needed """
		logging.debug('Monitoring constraints and params ')

	def _postprocessor(self) : 
		""" Run post processing steps """
		logging.debug("Running AR post processor")
		pass


	def solve(self) : 
		""" Runs the core routine """
		logging.debug('Beginning AR solver')	
		
		# Execute Preprocessing Steps
		self._preprocessor()
		
		# Execute the AR core
		self._core()

		# Run post processing steps
		self._postprocess()
	

	def _core(self) : 
		""" Core routine of the ARSolver class"""
		logging.debug('Beginning AR core')

		

class TrajectoryTracker(object) : 
	""" Tracks trajectory of the solution path """
	def __init__(self) : 
		pass

def _exp_splayer(ss) : 
	""" Splays the experiments and packs them into a new ss"""
	exp_list = ss.experiments
	for exp in exp_list : 
		ss_copy = copy.copy(ss)  # create a copy of original ss
		ss_copy.experiments = [exp]
		yield ss_copy



if __name__ == '__main__' :  
	pman = ParserManager()
	for ii,ss in enumerate(pman.get_gen_chou2006()) :
		# Run Alternating Regression and the result
		# Extract each individual ss and execute it 
		for expid,ss_exp in enumerate(_exp_splayer(ss)) : 
			print("Running ss: %s mod: %d exp: %d"% 
				(ss.name,ii,expid))	
			ar = ARSolver(ss_exp) 
			result_exp = ar.solve()


