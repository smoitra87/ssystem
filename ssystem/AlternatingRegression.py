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
import logging, copy, re, sys
from modifiers import ModifierChou2006
import numpy as np
import scipy as sp
import pylab as pl
import pdb


class ARSolver(object) :
	"""
Houses the Alternating Regression core routine and supporting scripts
to enforce constraints and good behavior of algorithm
	"""
	def __init__(self,ss) :  
		""" Init the AR solver""" 
		self.ss = ss
		self.logger = logging.getLogger('ss.ar')

	def _parse_initbound(self) : 
		""" Parse the soft constraints """
		logging.debug("Parsing initbound soft constraints")
		
	def _parse_initsol(self) : 
		""" Parse the initial solution """
		logging.debug("Parsing initsol initial solution")

		# Init initsol as an empty dict
		self.initsol = {}

		for varname in ['alpha','beta','g','h'] : 
			self._parse_var_initsol(varname)
	
	def _parse_var_initsol(self,varname) : 
		""" Apply default or specific values to attr with varname"""
		initsol = self.ss.constraint.initsol
		params = getattr(initsol,varname)
		nvars = len(self.ss.variables) # num of variables

		if varname in ('alpha','beta') : 
			self.initsol[varname] = np.ones(nvars)
			keys = params.keys()
			self.initsol[varname][:] = params['defaultInitialValue']
			for key in keys : 
				if re.match(varname+'_\d+',key)	:
					idx = int(key.split('_')[1])
					self.initsol[varname][idx-1] = params[key]
		elif varname in ('g','h') :
			self.initsol[varname] = np.ones([nvars,nvars])
			keys = params.keys()
			self.initsol[varname][:] = params['defaultInitialValue']
			for key in keys : 
				if re.match(varname+'_\d+_\d+',key)	:
					idr,idc = map(int,(key.split('_')[1:3]))
					self.initsol[varname][idr-1][idc-1] = params[key]
		
		else :
			logging.error("Unrecognized varname %s quitting.." \
			%(varname))
			sys.exit(1)


	def _parse_modelspace(self) :
		""" Parse the modelspace which are the hard constraints """
		logging.debug("Parsing modelspace hard constraints")	
	
		self.modelspace = {}
		
		for varname in ['alpha','beta','g','h'] : 
			self._parse_var_modelspace(varname)
	
	def _parse_var_modelspace(self,varname) : 
		""" Apply default or specific values to attr with varname"""

		modelspace = self.ss.constraint.modelspace
		params = getattr(modelspace,varname)
		nvars = len(self.ss.variables) # num of variables

		if varname in ('alpha','beta') : 
			keys = params.keys()
			var_range =  (params['defaultLowerBound'],\
				params['defaultUpperBound'])
			self.modelspace[varname] = [var_range]*nvars
			for key in keys : 
				if re.match(varname+'_\d+',key)	:
					idx = int(key.split('_')[1])				
					self.modelspace[varname][idx-1] = params[key]

		elif varname in ('g','h') :
			keys = params.keys()
			var_range = (params['defaultLowerBound'],\
				params['defaultUpperBound'])

			# This step is purely there cuz [[var_range]*nvars]*nvars
			# does not work
			varlist = []
			for ii in range(nvars) : 
				varlist.append([var_range]*nvars)
			self.modelspace[varname] = varlist
			for key in keys : 
				if re.match(varname+'_\d+_\d+',key)	:
					idr,idc = map(int,(key.split('_')[1:3]))
					self.modelspace[varname][idr-1][idc-1] = params[key]
		
		else :
			logging.error("Unrecognized varname %s quitting.." \
			%(varname))
			sys.exit(1)


	def _parse_softspace(self) :
		""" Compute the soft constraint space 
			Current strategy is to set it to 50% of hard constraint 
			limits
		"""
		self.softspace = {}
		
		for varname in ['alpha','beta','g','h'] : 
			self._parse_softspace_var(varname)
	
	def _parse_softspace_var(self,varname) : 
		""" Sets the values of alpha, beta, g and h from modelspace
			values
		"""
		modelspace = self.modelspace

		if varname in ('alpha','beta') : 
			vspace = modelspace[varname]
			vlist = []
			for v in vspace : 				
				if not type(v) is tuple : 
					vlist.append(v)
				else : 
					vlist.append((v[0]*0.5,v[1]*0.5))
			self.softspace[varname] = vlist
				
		elif varname in ('g','h') :
			vspace = modelspace[varname]
			vlist = []
			for vrow in vspace : 
				vlist2 = []
				for v in vrow : 				
					if not type(v) is tuple : 
						vlist2.append(v)
					else : 
						vlist2.append((v[0]*0.5,v[1]*0.5))
				vlist.append(vlist2)
			self.softspace[varname] = vlist
						
		else :
			logging.error("Unrecognized varname %s quitting.." \
			%(varname))
			sys.exit(1)

		
		

	def	_preprocessor(self) :  
		""" Run preprocessing on ss to make compatible with exp type
		"""
		logging.debug('Beginning preprocessor')
		
		# Parse entries from ss class
		self._parse_initsol()
		self._parse_modelspace()
		self._parse_softspace()
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


	def solve(self,maxiter=1000,tol=-7) : 
		""" Runs the core routine """
		logging.debug('Beginning AR solver')	
		
		# Execute Preprocessing Steps
		self._preprocessor()
		
		# Execute the AR core
		self._core()

		# Run post processing steps
		self._postprocessor()
	

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
