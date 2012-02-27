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


class ALRSolver(AR.ARSolver) :
	""" The ALR class is a variant of AR"""
	def __init__(self,_ss_dict) : 
		super(ALRSolver,self).__init__(_ss_dict)

	def _lasso(self) : 
		""" Makes the lasso call """
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
			alr = ALRSolver(ss_exp) 
			#result_exp = ar.solve(maxiter=1000,tol=10e-6)
			result_exp = alr.solve()

	
