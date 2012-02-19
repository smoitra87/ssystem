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
from utility import basedir,logdir,within_tuple
import logging, copy, re, sys, random
from modifiers import ModifierChou2006

import utility as util

import numpy as np
from numpy import linalg as LA
import scipy as sp
import pylab as pl

import pdb
from nose.tools import eq_
from other_models import Chou2006


class ARSolver(object) :
	"""
Houses the Alternating Regression core routine and supporting scripts
to enforce constraints and good behavior of algorithm
	"""
	def __init__(self,ss) :  
		""" Init the AR solver""" 
		self.ss = ss
		self.logger = logging.getLogger('ss.ar')
		
		# Run preprocessing steps
		self._preprocessor()
		

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

The following are the constraint concepts :-

Modelspace : Basically whatever is given in the S-system.
hardspace : A property of the algorithm which is used for 
catching crazy growing values
softspace  : Initially set to half modelspace and then adjusted to go
 up to modelspace

fullinfo:
    Regressors are not included
    Modelspace is enforced
    modelspace may or may not be suitably modified to encode up/down 
	regulation (+ve,-ve)
    Lenient modelspace is enforced as well

partialinfo:
    Union of prod degrad regressors
    Modelspace is enforced
    modelspace may or may not be suitably modified to encode up/down
	 regulation (+ve,-ve)
    Lenient modelspace enforced as well

noinfo:
    All regressors
    Modelspace is enforced
    modelspace may or may not be suitably modified to encode up/down 
	regulation (+ve,-ve)
    Lenient modelspace enforced as well
    
structure:
    All regressors
    Lenient modelspace enforced 

		"""
		logging.debug('Beginning preprocessor')
		
		# Parse entries from ss class
		self._parse_initsol()
		self._parse_modelspace()
		self._parse_softspace()
		self._parse_initbound()
		
		# Set regressors according to exptype
		self._set_regressors()

		# Deal with exptype ??
		if self.ss.exptype == 'fullinfo' :
			pass
		elif self.ss.exptype == 'partialinfo' : 
			pass
		elif self.ss.exptype == 'noinfo' : 
			pass
		elif self.ss.exptype == 'structure' :
			pass
		else : 
			logging.error('Did not recognize exptype %s'% \
			(self.ss.exptype))
			sys.exit(1)

		# Deal with equations
		self.equations = self.ss.equations

		# Set hard constraints - Constraints to catch crazily growing 
		# values, set by algorithm
		self.hardspace = {
		'alpha' : (0,100),
		'beta' : (0,100),
		'g' : (-10,-10),
		'h' : (-10,-10)
		}

		# Deal with noisy data ??

	def _set_regressors(self) : 
		""" Set the regressors according to exptype """
		logging.debug("Setting regressors")
	
		nvars = len(self.ss.variables)		

		if self.ss.exptype in ('noinfo','structure') : 
			prod = range(1,nvars+1)
			degrad = prod
			self.regressors = [{
				'prod' : prod,
				'degrad' : degrad
			}]*len(self.ss.equations) 

		elif self.ss.exptype == 'fullinfo' :
			self.regressors = []
			for eqn in self.ss.equations : 
				prod=  self._find_regressors(eqn,'g')
				degrad = self._find_regressors(eqn,'h')
				self.regressors.append({
					'prod':prod,
					'degrad':degrad
				})

		elif self.ss.exptype == 'partialinfo' : 
			self.regressors = []
			for eqn in self.ss.equations : 
				prod=  self._find_regressors(eqn,'g')
				degrad = self._find_regressors(eqn,'h')
				prod = list(set(prod).union(degrad))
				degrad = prod
				self.regressors.append({
					'prod':prod,
					'degrad':degrad
				})

		else :
			logging.error('Did not recognize exptype %s'% \
			(self.ss.exptype))
			sys.exit(1)

	def _find_regressors(self,eqn,varname) : 
		""" Find true regressors from eqn and variable"""
		true_params = self.ss._ss_params
		params = true_params[varname][eqn-1]		
		regressors = [ii+1 for ii,p in enumerate(params) if p!=0 ]		
		return regressors

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
				
		# Execute the AR core
		for exp in ss.experiments : 
			self._core(exp,maxiter,tol)

		# Run post processing steps
		self._postprocessor()
	

	def _core(self,exp,maxiter,tol) : 
		"""
Core routine of the ARSolver class

Note bd_i is [log(beta_i) hi1 .. hip]

		"""
		logging.debug('Beginning AR core')

		# Get profile from experiment
		prof = exp.profile	
		lp_list,ld_list,cp_list,cd_list = self._core_calc_design(prof)


		# Set initial solution params
		a_list,b_list,g_list,h_list = self._core_init_params()

		# Set Loop params
		self._continueLoop = True
		loopiter = 0
		
		for eqnid,eqn in enumerate(self.equations) :
			Cp = cp_list[eqnid]
			Cd = cd_list[eqnid]	
			Lp = lp_list[eqnid]			
			Ld = ld_list[eqnid]
			b = b_list[eqnid]
			h = h_list[eqnid]
			g = g_list[eqnid]
			a = a_list[eqnid]
			bd = np.concatenate(([np.log(b)],h))
			bp = np.concatenate(([np.log(a)],g))
			# Note this directly references slopes
			slopes = exp.profile.slopes[:,eqn-1]
			

			while(self._continueLoop == True) :
				
				# phase 1 components

				retcode,bp = self._core_phase1(slopes,Cp,bd,Ld,eqn)

				if retcode == 2 : 
					self.continueLoop = False
					self.logging.error(
					'Terminating eqn%d because of complex pain'%(eqn))
					break
				
				retcode = self._core_monitor_phase1(bp,eqnid,eqn)	
#	
#				# phase 2 components
#				a = a_list(eqnid-1)
#				g = g_list(eqnid-1)
#				bp = np.concatenate((np.log(a),g))
#				
#				retcode = self._core_phase2(Cd,slopes,bp)
#				
#				if retcode == BAD :
#					self.continueLoop = False
#					break
#	
#				self._core_monitor_phase1(bp,eqnid)			
#				
#				self.check_convergence()
				
				# For debug only
				self._continueLoop = False
	
	def _core_monitor_phase1(self,bp,eqnid,eqn) : 
		""" 
Monitor the bp params returned after phase1
Check that the bp params lie within the modelspace
if they don't fix them 
		"""
		reg_p = self.regressors[eqnid]['prod']
		a_mspace = self.modelspace['alpha'][eqn-1]
		g_mspace = self.modelspace['g'][eqn-1]
		g_mspace = [g_mspace[reg-1] for reg in reg_p]

		alpha = np.exp(bp[0])
		retcode = 0

		# check alpha_mspace
		if a_mspace == "nonZero"  :
			# Check if alpha is really close to zero
			if abs(np.exp(alpha)) < util.eps  :
				bp[0] = np.log(util.eps)
				retcode = 1
		else : 
			if not within_tuple(a_mspace,np.exp(alpha)) : 
				self.logger.debug(
				'alpha=%r out of mspace=%r'%(alpha,a_mspace))
				retcode = 1
				if alpha > a_mspace[1] : 
					bp[0] = np.log(a_mspace[1]-util.eps)
				else : 
					bp[0] = np.log(a_mspace[0]+util.eps)
		
		# check g_mspace. Loop through all the g's there
		#	pdb.set_trace()
		for ii,g in enumerate(bp[1:]) : 
			g_space = g_mspace[ii]
			if g_space == "nonZero"  :
				# Check if param is really close to zero
				if abs(g) < util.eps  :
					bp[1+ii] = random.choice([-util.eps,util.eps])
					retcode = 1
			else : 
				if not within_tuple(g_space,g) : 
					self.logger.debug(
					'g%d=%r out of mspace=%r. Fixing..'%
					(ii+1,g,g_space))
					retcode = 1
					if alpha > a_mspace[1] : 
						bp[ii+1] = g_space[1]-util.eps
					else : 
						bp[ii+1] = g_space[0]+util.eps
				
		
		return retcode				


	
	def _core_phase1(self,slopes,Cp,bd,Ld,eqn)  : 
		""" 
Run phase1  : 
	Calculates degrad terms and estimates bp
		Return codes:
		0 : Success
		1 : Complex trouble, but solved
		2 : Complex trouble, but could not solve
		"""

		degrad = self._core_calc_degrad(bd,Ld)
		yd_ = slopes + degrad
		retcode = 0 # Assuming success

		# If yd_ <= 0 . Try to solve it by fixing beta 
		while (yd_ <= 0.0).any() and \
		np.exp(bd[0])< ar.modelspace['beta'][eqn-1][1] :

			self.logger.debug(\
			"Found complex pain. Dealing with it..")
			retcode =  1 # Found complex pain
			
			bd[0] += np.log(2) # Mult beta by 2 to fix complex pain
			degrad = self._core_calc_degrad(bd,Ld)
			yd_ = slopes+degrad

		if (yd_ <= 0.0).any() : 
			self.logger.error('Could not deal with complex pain')
			retcode = 2
			return retcode,None	
					
		
		yd = np.log(yd_) 
		bp = np.dot(Cp,yd)

		return retcode,bp	
	


	def _core_calc_degrad(self,bd,Ld) : 
		""" 
Calculate the degradation values
Given bd = [log(b) hi1 hi2 .. hip]
Ld = Design matrix
It returns [n_sample x 1] array of degrad values
	
		"""
		degrad = np.dot(Ld,bd) # Do matrix multiplication 
		degrad = np.exp(degrad) # Exponentiate to convert log to real
		return degrad

	def _core_calc_prod(self,bp,Lp) : 
		"""
Calculate the production values
Given bp = [log(a) gi1 gi2 .. gip]
Lp = Design matrix
It returns [n_sample x 1] array of prod values
		"""
		prod = np.dot(Lp,bp)
		prod = np.exp(prod)
		return prod


	def _core_phase2(self,Cd,slopes) :
		""" Compute phase 2 """
		pass

	def _core_monitor_phase2(self,bd) : 
		""" Monitor the bd params returned after phase2 """
		pass

	def _core_init_params(self) : 
		""" 
Initialize the params for b and h from initsol, for all other params 
set some dummy values. Separate every param equation wise and stick it 
in a list
		"""
		a_list,b_list = [],[]
		g_list,h_list = [],[]
		
		
		for eqnid,eqn in enumerate(self.equations) : 
			reg_p = self.regressors[eqnid]['prod']
			reg_d = self.regressors[eqnid]['degrad']
			h_eqn = self.initsol['h'][eqn-1]

			a_list.append(1.0) # dummy
			b_list.append(self.initsol['beta'][eqn-1])
		
			g_list.append(np.zeros(len(reg_p)))
		
			h_eqn = np.array([h_eqn[reg-1] for reg in reg_d])
			h_list.append(h_eqn)
	
		return (a_list,b_list,g_list,h_list)
	
	def _core_calc_design(self,prof) : 
		""" 
Calculate design matrices and return list	
Calculate Lp, Ld, Cp and Cd 

lp_list contains a list of Lp
Lp= |1 log(X1(t1)) . . log(Xp(t1)) |
	|1 log(X1(t2))				   |
	|.		.					   |
	|.		.					   |
	|1 log(X1(tn)) . . log(Xp(tn)) |

Here 1..p may depend on what regressors are actually being selected
according to exptype

Cp = inv(Lp'*Lp)*Lp'

Similarly for Ld

		 """
		lp_list,ld_list = [],[]
		cp_list,cd_list = [],[]
		
		
		for eqnid,eqn in enumerate(self.equations) : 
			reg_p = self.regressors[eqnid]['prod']
			reg_d = self.regressors[eqnid]['degrad']
	
			Lp = np.ones(prof.n_sample)
			Ld = np.ones(prof.n_sample)
		
			# Get regressor values
			X_p = [np.log(prof.var[:,reg-1]) for reg in reg_p ]
			X_d = [np.log(prof.var[:,reg-1]) for reg in reg_d ]
			
			Lp = np.vstack((Lp,np.array(X_p))).T
			Ld = np.vstack((Ld,np.array(X_d))).T

			# Calculate Cp
			Cp = np.dot(LA.inv(np.dot(Lp.T,Lp)),Lp.T)
			Cd = np.dot(LA.inv(np.dot(Ld.T,Ld)),Ld.T)


			# Append Lp,Ld,Cp and Cd to relevant lists
			lp_list.append(Lp)
			ld_list.append(Ld)
			cp_list.append(Cp)
			cd_list.append(Cd)

			
		return (lp_list,ld_list,cp_list,cd_list)


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


	
