""" 
This module contains several classes implementing a number of Feature
Selection methods for S-systems. Currently the purpose is to service
AR like algorithms. 

The host of classes implemented here are : 
FeatSel - Basic skeleton class that needs to be derived in order to 
build a concrete Feature Selection class
LassoARFeatSel - Does Feature Selection by running Lasso in an AR like
scheme

author: Subhodeep Moitra(smoitra@cs.cmu.edu)

BSD License..!
"""

import sys, os, pdb, logging

import numpy as np
import scipy as sci
import pylab as pl
from numpy import linalg as LA 


import utility as util
from utility import overrides

class FeatSel(object) : 
	""" Init the Feat Sel class"""
	def __init__(self) : 
		self.features = None
		self.params = None

	def _preprocessor(self) : 
		""" Run preprocessor for Feature Selection, 
			This is a virtual class. Make it concrete..!
		"""
		pass

	def find_features(self) : 
		""" Find the relevant features
			This is a virtual class. Make it concrete..!
		"""
		pass


class LassoARFeatSel(FeatSel) : 
	""" Perform AR FeatSel by doing Lasso"""
	def __init__(self,ar,expid,**kwargs) :
		self.logger = logging.getLogger("ss.ar.fsl")
		self.name = "FSARL"  
		self.ar = ar
		self.expid = expid
		self.exp = self.ar.ss.experiments[expid]
		self.prof = self.exp.profile
		self.paramspace = None		
		self.regressors = None
		self.features = None
		self.regfunc = self._least_squares
		self._preprocessor(**kwargs)

	@overrides(FeatSel)
	def _preprocessor(self,**kwargs) : 
		"""	
		Run preprocessor for Lasso AR Feature Selection 
		"""
		# Parse the keyword args that core needs
		self._core_ar_kwarg(**kwargs)

		# Get profile from experiment
		self.lp_list,self.ld_list,self.cp_list,self.cd_list = \
			self._core_calc_design(self.prof)

	def _Iregfunc_handler(self,L,C,y) : 
		""" Selects whether to send L,C depending on AR/ALR 
			This acts like an interface that is defined by AR/ALR and 
			other variants
		"""
		if self.name == "FSARL" : 
			return self.regfunc(C,y)
		else :
			self.logger.error("Unknown name %r"%(self.name))
			sys.exit(1)

	def _least_squares(self,C,y) :
		""" Simple linear regression to calc bp """
		b = np.dot(C,y)
		return b	

	def gen_rand_modelspace(self) : 
		""" Returns random beta and h values from modelspace """
		rand_param = {'beta':[],'h':[]}
		msp = self.ar.modelspace
		for eqnid, eqn in enumerate(self.ar.equations) : 
			m_beta = msp['beta'][eqnid]
			low,high = m_beta[0],m_beta[1]
			rand_param['beta'].append(\
				np.random.triangular(low,10,high))
			m_h = msp['h'][eqnid]
			rand_param['h'].append([])
			for m_h_var in m_h : 
				low,high = m_h_var[0],m_h_var[1]
				rand_param['h'][eqnid].append(\
					np.random.triangular(low,0,high))
		return rand_param

	@overrides(FeatSel)
	def find_features(self,**kwargs) : 
		""" Find the relevant features using Lasso AR """
		self.features = []
		for eqnid,eqn in enumerate(self.ar.equations) :

			feats,rc1,rc2 = self.find_feat_eqn(eqnid,eqn,**kwargs)
			feats = {'prod':[3],'degrad':[1,2]}
			self.features.append(feats)
	
		return self.features

	def find_feat_eqn(self,eqnid,eqn,**kwargs) : 
		"""
		Find the features for a given equation
		"""
		Cp = self.cp_list[eqnid]
		Cd = self.cd_list[eqnid]	
		Lp = self.lp_list[eqnid]			
		Ld = self.ld_list[eqnid]
		slopes = self.prof.slopes[:,eqn-1]
		nScan = 100 # num times to seed rand initsol

		for scanid in xrange(nScan) : 
			self.art = ARTracker(ar=self,eqn=eqn,maxiter=1,**kwargs)
			param = self.gen_rand_modelspace()
			b = param['beta'][eqnid]
			h = np.array(param['h'][eqnid])
			a = b # dummy
			g = h # dummy
			bd = np.concatenate(([np.log(b)],h))
			bp = np.concatenate(([np.log(a)],g))
			# Note this directly references slopes
			retcode1 = retcode2 = 0	
	
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
	
				#------------------------------------------------------
				# PHASE II
				#------------------------------------------------------
				
				retcode2,arp.bd,ssed = self._core_phase2(arp,eqn)
	
				if retcode2 == 2 : 
					self.art.continueLoop = False
					self.logger.warning(
					'Terminating eqn%d because of complex pain'%(eqn))
					break
	
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
	
			if retcode1 == 2 or retcode2 == 2 :
				continue
			


	def _core_ar_kwarg(self,**kwargs) : 
		""" Parses kwargs for core module in ARSolver"""
		pass
	
	
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
		
		nvars = len(self.ar.ss.variables)
		
		for eqnid,eqn in enumerate(self.ar.equations) : 
			reg_p = range(1,nvars+1)
			reg_d = range(1,nvars+1)
			
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

		if (yd_ <= 0.0).any()  : 
			retcode = 2
			return retcode,None,None
							
		yd = np.log(yd_) 
		bp = self._Iregfunc_handler(Lp,Cp,yd)
		#bp = np.dot(Cp,yd)
	
		# Calculate ssep
		ssep = self._core_calc_sse(yd,Lp,bp)

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

		if (yp_ <= 0.0).any() :
			retcode = 2
			return retcode,None,None
	
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
		


if __name__ == '__main__' : 
	fsel = LassoARFeatSel()
	fsel.find_features()


	
