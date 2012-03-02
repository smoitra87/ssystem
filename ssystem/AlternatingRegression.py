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
from FeatureSelector import LassoARFeatSel


import utility as util

import numpy as np
from numpy import linalg as LA
import scipy as sp
import pylab as pl

import pdb
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
		self.regfunc = self._least_squares
		self.name = "AR"	
		self.exptype = ss.exptype	
		self._IARTracker = ARTracker
		self.regressors = None
		self.regressors_true = None
		self.modelspace = None
		self.art = None
		self.exp_art = None
		self.all_exp_art = None
		self.equations = None
		self._regressors_fix = None
		self._core_Ikwarg = self._core_ar_kwarg

		# Run preprocessing steps
		self._preprocessor()
			
	def _Iregfunc_handler(self,L,C,y) : 
		""" Selects whether to send L,C depending on AR/ALR 
			This acts like an interface that is defined by AR/ALR and 
			other variants
		"""
		if self.name == "AR" : 
			return self.regfunc(C,y)
		else :
			self.logger.error("Unknown name %r"%(self.name))
			sys.exit(1)

	def _least_squares(self,C,y) :
		""" Simple linear regression to calc bp """
		b = np.dot(C,y)
		return b	

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

	def	_preprocessor(self) :  
		""" Run preprocessing on ss to make compatible with exp type

The following are the constraint concepts :-

Modelspace : Basically whatever is given in the S-system.

fullinfo:
    Regressors are not included

partialinfo:
    Union of prod degrad regressors

noinfo:
    All regressors
   
structure: --- Not implemented yet
    All regressors
    Lenient modelspace enforced 

		"""
		logging.debug('Beginning preprocessor')
		
		# Parse entries from ss class
		self._parse_initsol()
		self._parse_modelspace()
		self._parse_initbound()
		
		# Set regressors according to exptype
		self._set_regressors()

		# Deal with equations
		self.equations = self.ss.equations

		# Deal with noisy data ??

	def _set_regressors(self) : 
		""" Set the regressors according to exptype """
		logging.debug("Setting regressors")
	
		nvars = len(self.ss.variables)		

		self.regressors_true = []
		for eqn in self.ss.equations : 
			prod=  self._find_regressors(eqn,'g')
			degrad = self._find_regressors(eqn,'h')
			self.regressors_true.append({
				'prod':prod,
				'degrad':degrad
			})

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

		# The regressors that are to be artificially fixed at every 
		# iteration
		# This is only necessary for the noinfo and partial info case
		self._regressors_fix = []
		for eqn in self.ss.equations : 
			prod=  self._find_regressors(eqn,'g')
			degrad = self._find_regressors(eqn,'h')
			prod_all = degrad_all = range(1,nvars+1)
			prod = list(set(prod_all).difference(prod))
			degrad = list(set(degrad_all).difference(degrad))
			self._regressors_fix.append({
				'prod':prod,
				'degrad':degrad
			})

	def _find_regressors(self,eqn,varname) : 
		""" Find true regressors from eqn and variable"""
		true_params = self.ss._ss_params
		params = true_params[varname][eqn-1]		
		regressors = [ii+1 for ii,p in enumerate(params) if p!=0 ]		
		return regressors

	def _core(self,exp,**kwargs) : 
		"""
Core routine of the ARSolver class
Note bd_i is [log(beta_i) hi1 .. hip]
		"""
		logging.debug('Beginning AR core')
	
		# Parse the keyword args that core needs
		self._core_Ikwarg(**kwargs)

		# Get profile from experiment
		prof = exp.profile	
		lp_list,ld_list,cp_list,cd_list = self._core_calc_design(prof)

		# Set initial solution params
		a_list,b_list,g_list,h_list = self._core_init_params()

		# Set Loop params
		
		for eqnid,eqn in enumerate(self.equations) :
			self.art = self._IARTracker(ar=self,eqn=eqn,**kwargs)
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
			
			arp = ARParams(Cp,Cd,Lp,Ld,b,h,g,a,bd,bp,slopes)

			while(self.art.continueLoop == True) :				
				#------------------------------------------------------
				# PHASE I
				#-----------------------------------------------------

				retcode1,arp.bp,ssep = self._core_phase1(arp,eqn)

				if retcode1 == 2 : 
					self.art.continueLoop = False
					self.logger.warning(
					'Terminating eqn%d because of complex pain'%(eqn))
					break			

				# fix params after phase1	
				arp.bp = self._core_fix_params(arp.bp,\
					phase=1,eqnid=eqnid)
				
				# Monitor and fix bp params if needed	
				#retcode1 = self._core_monitor(bp,eqnid,eqn,phase=1)	
				#------------------------------------------------------
				# PHASE II
				#------------------------------------------------------
				
				retcode2,arp.bd,ssed = self._core_phase2(arp,eqn)

				if retcode2 == 2 : 
					self.art.continueLoop = False
					self.logger.warning(
					'Terminating eqn%d because of complex pain'%(eqn))
					break

				# fix params after phase2
				arp.bd = self._core_fix_params(arp.bd,\
					phase=2,eqnid=eqnid)

				# Monitor bd , may be used at a later stage 
				#retcode2 = self._core_monitor(bd,eqnid,eqn,phase=2)					
				# Calculate overall sse for fit to curve
				sse = self._core_calc_sse2(arp)		
	
				# Perform bookkeeping and convergence checks
				self.art.set_SSE(ssep,ssed,sse)
				self.art.set_params(arp.bp,arp.bd)
				self.art.set_retcodes(retcode1,retcode2)
				self.art.bookkeep()	
				self.art.check_termination()			
	
			# store the alpha and the beta values
			# Append AR Trace to Experiment Trace list
			bp,bd = arp.bp,arp.bd
			if not self.art.save_trace : 
				if (bp is None) or (bd is None): 
					_params = None	
				else : 
					_params = {}
					_params['alpha'] = np.exp(bp[0])
					_params['beta'] = np.exp(bd[0])
					_params['g'] = bp[1:]
					_params['h'] = bd[1:]
				self.art.params = [_params]
			
			self.exp_art['eqns'].append(self.art)

			# Don't store results if complex pain could not be solved
			if retcode1 == 2 or retcode2 == 2 : 
				self.logger.info("Premature Failure ! Skipping equation						..")
				continue

			# Logging information after a particular equation finishes
			if(self.art.converged == True ) :
				self.logger.info("Convergence Succeeded..!")
				self.logger.debug("ssep=%f,ssed=%f,sse=%f"%\
					(ssep,ssed,sse))
				self.logger.debug("alpha=%r"%(np.exp(bp[0])))
				self.logger.debug("g=%r"%(bp[1:]))
				self.logger.debug("beta=%r"%(np.exp(bd[0])))
				self.logger.debug("h=%r"%(bd[1:]))
				#util.plot_pair(slopes,prod-degrad,['true','estimate'])
			else : 
				self.logger.info("Convergence Failed..!")
				self.logger.debug("maxiter_exceeded ? %r"%\
					(self.art.maxiter_exceeded))
				self.logger.debug("ssep=%f,ssed=%f,sse=%f"%\
					(ssep,ssed,sse))
				self.logger.debug("alpha=%r"%(np.exp(bp[0])))
				self.logger.debug("g=%r"%(bp[1:]))
				self.logger.debug("beta=%r"%(np.exp(bd[0])))
				self.logger.debug("h=%r"%(bd[1:]))
				#util.plot_pair(slopes,prod-degrad,['true','estimate'])

	def _core_ar_kwarg(self,**kwargs) : 
		""" Parses kwargs for core module in ARSolver"""
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
			g_eqn = self.initsol['g'][eqn-1]


			a_list.append(self.initsol['alpha'][eqn-1])
			b_list.append(self.initsol['beta'][eqn-1])
			
			g_eqn = np.array([g_eqn[reg-1] for reg in reg_p])
			h_eqn = np.array([h_eqn[reg-1] for reg in reg_d])
			h_list.append(h_eqn)
			g_list.append(g_eqn)
	
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

	def _core_fix_params(self,bx,phase,eqnid) : 
		""" Fixes params to zero depending on nature of problem """
		if self.ss.exptype == 'fullinfo' :
			return bx
		elif self.ss.exptype in ('noinfo','partialinfo') :
			reg_fix = self._regressors_fix[eqnid] 
			if phase == 1 : 
				prod_fix = reg_fix['prod']
				for varid in prod_fix : 
					bx[varid] = 0.0
			elif phase == 2 : 			
				degrad_fix = reg_fix['degrad']
				for varid in degrad_fix :
					bx[varid] = 0.0 
			else : 
				self.logger.error("Incorrect phase %r"%(phase))
				sys.exit(1)
		else :
			self.logger.error("Unrecognized exptype %s quitting..."%\
			(self.ss.exptype))

		return bx

	def _core_phase1(self,arp,eqn)  : 
		""" 
Run phase1  : 
	Calculates degrad terms and estimates bp
		Return codes:
		0 : Success
		1 : Complex trouble, but solved
		2 : Complex trouble, but could not solve
		"""
		Cp = arp.Cp
		slopes = arp.slopes
		Lp  = arp.Lp
		Ld = arp.Ld
		bd = arp.bd

		degrad = self._core_calc_degrad(bd,Ld)
		yd_ = slopes + degrad
		retcode = 0 # Assuming success

		# If yd_ <= 0 . Try to solve it by fixing beta 
		while (yd_ <= 0.0).any() and \
		2*np.exp(bd[0])< self.modelspace['beta'][eqn-1][1] :

			self.logger.debug(\
			"Found complex pain in phase 1. Dealing with it..")
			retcode =  1 # Found complex pain
			
			bd[0] += np.log(2) # Mult beta by 2 to fix complex pain
			degrad = self._core_calc_degrad(bd,Ld)
			yd_ = slopes+degrad

		if retcode == 1: 
			self.logger.debug("iter=%d,fixed beta=%f"%\
			(self.art.loopiter,np.exp(bd[0])))

		if (yd_ <= 0.0).any() : 
			self.logger.warning(\
			'Could not deal with complex pain in phase1')
			self.logger.debug("iter=%d,beta=%f"%\
				(self.art.loopiter,np.exp(bd[0])))

			retcode = 2
			return retcode,None,None
							
		yd = np.log(yd_) 
		bp = self._Iregfunc_handler(Lp,Cp,yd)
		#bp = np.dot(Cp,yd)
	
		# Calculate ssep
		ssep = self._core_calc_sse(yd,Lp,bp)
		# DEBUG
		#util.plot_pair(yd,np.dot(Lp,bp),labels=['yd','Lp*bp'])			
		return retcode,bp,ssep	

	def _core_phase2(self,arp,eqn)  : 
		""" 
Run phase2  : 
	Calculates prod terms and estimates bd
		Return codes:
		0 : Success
		1 : Complex trouble, but solved
		2 : Complex trouble, but could not solve
		"""
		slopes = arp.slopes
		Cd = arp.Cd
		Lp = arp.Lp
		Ld = arp.Ld
		bp = arp.bp

		prod = self._core_calc_prod(bp,Lp)
		yp_ =  prod - slopes
		retcode = 0 # Assuming success

		# If yp_ <= 0 . Try to solve it by fixing alpha 
		while (yp_ <= 0.0).any() and \
		2*np.exp(bp[0])< self.modelspace['alpha'][eqn-1][1] :

			self.logger.debug(\
			"Found complex pain in phase2. Dealing with it..")
			retcode =  1 # Found complex pain
			
			bp[0] += np.log(2) # Mult alpha by 2 to fix complex pain
			prod = self._core_calc_prod(bp,Lp)
			yp_ = prod - slopes

		if (yp_ <= 0.0).any() : 
			self.logger.warning(\
				'Could not deal with complex pain in phase2')
			self.logger.debug("iter=%d,alpha=%f"%\
				(self.art.loopiter,np.exp(bp[0])))
			retcode = 2
			return retcode,None,None	
					
		if retcode == 1: 
			self.logger.debug("iter=%d,fixed alpha=%f"%\
			(self.art.loopiter,np.exp(bp[0])))
	
		yp = np.log(yp_) 
		bd = self._Iregfunc_handler(Ld,Cd,yp)
		#bd = np.dot(Cd,yp)

		# Calculate ssed
		ssed = self._core_calc_sse(yp,Ld,bd)
		#util.plot_pair(yp,np.dot(Ld,bd),labels=['yp','Ld*bd'])

		return retcode,bd,ssed


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

	def _core_calc_sse2(self,arp) : 
		""" Calculates SSE for overall fit to curve. 
			The other SSE function calculates for individual ssep/d
		"""
		prod = self._core_calc_prod(arp.bp,arp.Lp)
		degrad = self._core_calc_degrad(arp.bd,arp.Ld)
		e = arp.slopes - prod + degrad
		sse = np.dot(e,e)			
		return sse

	def _core_calc_sse(self,y,L,b)  : 
		""" Calculates SSE """
		e = y - np.dot(L,b)
		sse = np.dot(e,e)
		return sse
	
	def _postprocessor(self) : 
		""" Run post processing steps """
		logging.debug("nning AR post processor")
		pass

	def solve(self,**kwargs) : 
		""" Runs the core routine """
		logging.debug('Beginning AR solver')	
	
		# List of trackers for all experiments
		self.all_exp_art = []	

		# Execute the AR core
		for expid,exp in enumerate(self.ss.experiments) : 
				self.exp_art = {} # AR tracker for a single experiments
				self.exp_art['id'] = expid+1
				self.exp_art['eqns'] = []
				self._core(exp,**kwargs)
				self.all_exp_art.append(self.exp_art)

		#Run post processing steps
		self._postprocessor()

	#-----------------------------------------------------------------
	# Depreceated Code Section
	#-----------------------------------------------------------------

	# Depreceated...!!! :-(
	def _core_monitor(self,bx,eqnid,eqn,phase) : 
		""" 
Here bx stands for either bp or bd, depending on which phase called
Monitor the bx params returned after respective phases
Check that the bx params lie within the modelspace
if they don't fix them 
		"""
		# Depending on which phase is requesting monitor
		if phase == 1 :
			dynamics = 'prod'
			var1 = 'alpha'
			var2 = 'g'
		elif phase == 2 : 
			dynamics = 'degrad'
			var1 = 'beta'
			var2 = 'h'
		else : 
			self.logger.error("Unknown phase:%r"%(phase))
			sys.exit(1)

		reg_x = self.regressors[eqnid][dynamics]
		ab_mspace = self.modelspace[var1][eqn-1]
		gh_mspace = self.modelspace[var2][eqn-1]
		gh_mspace = [gh_mspace[reg-1] for reg in reg_x]

		ab = np.exp(bx[0])
		retcode = 0

#		print "reg_x is", reg_x
#		print "ab_mspace is", ab_mspace
#		print "gh_mspace is", gh_mspace
#		print "ab is", ab
#		print "bx is", bx

		# check abb_mspace
		if ab_mspace == "nonZero"  :
			# Check if ab is really close to zero
			if abs(np.exp(ab)) < util.eps  :
				bx[0] = np.log(util.eps)
				retcode = 1
		else : 
			if not within_tuple(ab_mspace,ab) : 
				self.logger.debug(
				'Monitor found violation after phase%d',phase)
				self.logger.debug(
				'%s=%r out of mspace=%r'%(var1,ab,ab_mspace))
				retcode = 1
				if ab > ab_mspace[1] : 
					bx[0] = np.log(ab_mspace[1]-util.eps)
				else : 
					bx[0] = np.log(ab_mspace[0]+util.eps)
		
		# check gh_mspace. Loop through all the g's there
		#	pdb.set_trace()
		for ii,gh in enumerate(bx[1:]) : 
			gh_space = gh_mspace[ii]
			if gh_space == "nonZero"  :
				# Check if param is really close to zero
				if abs(gh) < util.eps  :
					bx[1+ii] = random.choice([-util.eps,util.eps])
					retcode = 1
			else : 
				if not within_tuple(gh_space,gh) : 
					self.logger.debug(
				'Monitor found violation after phase%d',phase)
					self.logger.debug(
					'%s_%d=%r out of mspace=%r. Fixing..'%
					(var2,ii+1,gh,gh_space))
					retcode = 1
					if gh > gh_space[1] : 
						bx[ii+1] = gh_space[1]-util.eps
					else : 
						bx[ii+1] = gh_space[0]+util.eps
				
		
		return retcode				


class ARTracker(object) : 
	""" 
Contains values for tracking convergence, maxiter, tolerance
	"""
	def __init__(self,ar,eqn,save_trace=False,\
			maxiter=10000,tol=10e-6,**kwargs) : 
		self.maxiter = maxiter
		self.ar = ar
		self.eqn = eqn # The eqn which is being tracked
		self.tol = tol
		self.continueLoop = True
		self.loopiter = 0
		self.ssep = []
		self.ssed = []
		self.sse = []
		self.rc1 = []
		self.rc2 = []
		self.params = [] # Contains params from iterations
		self.converged = False
		self.maxiter_exceeded = False
		self.save_trace = save_trace
		self.sseMax = 15

	def bookkeep(self) : 
		""" Perform bookkeeping operations for AR algo"""
		self.loopiter += 1

	def check_termination(self) : 
		""" Check convergence for algo"""
		rc_all = self.rc1[-1]+self.rc2[-1]
		if rc_all == 0 and self.ssep[-1] < self.tol \
			 and self.ssed[-1] < self.tol : 
			self.converged = True
			self.continueLoop = False 
	
		if self.loopiter >= self.maxiter : 
			self.maxiter_exceeded = True
			self.converged = False
			self.continueLoop = False
	
	# Check any obvious SSE inconsistencies
		if self.ssep[-1] > self.sseMax or self.ssed[-1] > self.sseMax : 
			self.continueLoop = False							

	def set_SSE(self,ssep,ssed,sse) :
		""" Set the SSE values for the diff phases"""
		if self.save_trace :
			# save the ret codes from the last iter
			self.ssep.append(ssep)
			self.ssed.append(ssed)
			self.sse.append(sse)
		else : 
			self.ssep = [ssep]
			self.ssed = [ssed]
			self.sse = [sse]

	def set_retcodes(self,rc1,rc2) : 
		""" Set the return codes for phase I and II"""
		if self.save_trace : 
			self.rc1.append(rc1)
			self.rc2.append(rc2)
		else : 
			# Save the return codes from the last iter
			self.rc1 = [rc1]
			self.rc2 = [rc2]

	def set_params(self,bp,bd) : 
		""" Save the params """
		if self.save_trace : 
			if (bp is None) or (bd is None): 
				_params = None	
			else : 
				_params = {}
				_params['alpha'] = np.exp(bp[0])
				_params['beta'] = np.exp(bd[0])
				_params['g'] = bp[1:]
				_params['h'] = bd[1:]
	
			self.params.append(_params)

	def plot_params(self,sel='all',niter=None) : 
		""" Plot all the params 

		"""
		alpha = np.array([param['alpha'] for param in self.params])
		beta = np.array([param['beta'] for param in self.params])
		g = np.array([param['g'] for param in self.params])
		h = np.array([param['h'] for param in self.params])
		pdict = {
			'alpha' : alpha,
			'beta' : beta,
			'g' : g,
			'h' : h
		}	

		# Plot true params
		if sel == 'all' : 
			for key in pdict.keys() : 
				key_t = self.ar.ss._ss_params[key][self.eqn-1]
				if type(key_t) is list : 
					for t in key_t : 
						pl.axhline(y=t,color="black",linewidth=2,\
						linestyle="dashed")
				else : 
					pl.axhline(y=key_t,color="black",linewidth=2,\
						linestyle="dashed")
		else : 
			key_t = self.ar.ss._ss_params[sel][self.eqn-1]
			if type(key_t) is list : 
				for t in key_t : 
					pl.axhline(y=t,color="black",linewidth=2,\
						linestyle="dashed")
			else : 
				pl.axhline(y=key_t,color="black",linewidth=2,\
						linestyle="dashed")


		# Plot estimated params
		if sel == 'all' : 
			for key in pdict.keys() : 
				pl.plot(pdict[key],label=key)
			pl.legend()
			pl.xlabel("Iteration")
			pl.ylabel("Params")
		else : 
			pl.plot(pdict[sel],label=sel)
			pl.legend()
			pl.xlabel("Iteration")
			pl.ylabel(sel)
		pl.show()

	def plot_sse(self,sel='pd',niter=None,log=True) : 
		""" Plot the sse values """
		if not niter :
			niter = len(self.sse)

		if log : 
			sse = np.log(self.sse)
			ssep = np.log(self.ssep) 
			ssed = np.log(self.ssed)
		else : 
			sse = self.sse
			ssep = self.ssep
			ssed = self.ssed		

		if sel == "all" :
			pl.plot(range(niter),sse[:niter],label="SSE")
			pl.plot(range(niter),ssep[:niter],label="SSEp")
			pl.plot(range(niter),ssed[:niter],label="SSEd")
		else : 
			pl.plot(range(niter),ssep[:niter],label="SSEp")
			pl.plot(range(niter),ssed[:niter],label="SSEd")

		pl.xlim(-5,niter+5)
		if log : 
			pl.ylabel('log SSE')
		else : 
			pl.ylabel('SSE')
		pl.xlabel('niter')		

		pl.legend(loc="upper right")
		pl.show()

def _exp_splayer(ss) : 
	""" Splays the experiments and packs them into a new ss"""
	exp_list = ss.experiments
	for exp in exp_list : 
		ss_copy = copy.copy(ss)  # create a copy of original ss
		ss_copy.experiments = [exp]
		yield ss_copy

class ARParams(object) : 
	""" This object stores the AR params currently in use"""
	def __init__(self,Cp,Cd,Lp,Ld,b,h,g,a,bd,bp,slopes) : 
		self.Cp = Cp
		self.Cd = Cd
		self.Lp = Lp
		self.Ld = Ld
		self.b = b
		self.h = h
		self.g = g
		self.a = a
		self.bd = bd
		self.bp = bp
		self.slopes = slopes
		

class ARFeatLassoSolver(ARSolver) : 
	""" Lasso Feature Selector for doing AR"""
	def __init__(self,ss) :
		""" Init for Lasso Feat Solver""" 
		self.logger = logging.getLogger("ss.ar.fl")
		self.logger.debug("Calling Constructor for AR")
		super(ARFeatLassoSolver,self).__init__(ss)



if __name__ == '__main__' :  
	pman = ParserManager()
	for ii,ss in enumerate(pman.get_gen_chou2006()) :
		# Run Alternating Regression and the result
		# Extract each individual ss and execute it 
		for expid,ss_exp in enumerate(_exp_splayer(ss)) : 
			print("Running ss: %s mod: %d exp: %d"% 
				(ss.name,ii,expid))	
			ar = ARFeatLassoSolver(ss_exp)
			result_exp = ar.solve(save_trace=False,maxiter=10000,tol=10e-5)
			#result_exp = ar.solve()

	
