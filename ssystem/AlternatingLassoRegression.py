""" 
This module implements the Alternating Lasso Regression 

Input  : 
An S-system model encoded with the type of experiment encoded

@author: Subhodeep Moitra (smoitra@cs.cmu.edu)

BSD License

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
import AlternatingRegression as AR
from AlternatingRegression import ARSolver

from sklearn.linear_model import Lasso,LassoCV,LarsCV,LassoLars,\
	LassoLarsCV


class ALRSolver(AR.ARSolver) :
	""" The ALR class is a variant of AR"""
	def __init__(self,_ss_dict) : 
		super(ALRSolver,self).__init__(_ss_dict)
		self.name = "ALR"
		self.regfunc = self._lasso
		self.l1penalty = 0.0

	def _Iregfunc_handler(self,L,C,y) : 
		""" Selects whether to send L,C depending on AR/ALR """
		if self.name == "ALR" : 
			return self.regfunc(L,y)
		else :
			self.logger.error("Unknown name %r"%(self.name))
			sys.exit(1)

	def _lasso(self,X,y) : 
		""" Makes the lasso call 
			
		Parameters 
		----------
		X : The Lp or Ld matrix
		y : The yp or yd response vector
		
		Here alpha stands for the L1 penalty applied to lasso

		"""
		
		lasso = Lasso(alpha=self.l1penalty)
		lasso_model = lasso.fit(X,y)
		gh = lasso_model.coef_
		gh = gh[1:] # Throw out the first coeff since it is intercept
		ab = lasso_model.intercept_
		bx = np.concatenate(([ab],gh))
		return bx


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
			alr = ALRSolver(ss_exp) 
			result_exp = alr.solve(maxiter=10000,tol=10e-6)
			#result_exp = alr.solve()

	
