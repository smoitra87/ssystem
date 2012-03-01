"""
Test cases for Alternating Regression Module
"""

from nose.tools import *

from ssystem.AlternatingRegression import ARSolver
from ssystem.AlternatingLassoRegression import ALRSolver
from ssystem.other_models import Chou2006
import numpy.testing as npt
import numpy as np
from numpy import linalg as LA

from ssystem.utility import same_dict
from ssystem.base import Profile

from nose.plugins.attrib import attr


class TestALRSolverChou2006(object) :
	""" Test cases for the ALR solver """
	def setUp(self) : 
		self.ss = Chou2006()
		self.ar = ARSolver(self.ss)
	
	def tearDown(self) :
		pass


	@attr('slow')
	def test_fullinfo1(self) : 
		""" RUN ALR (slow) eqn 1 under full info and alpha=0.0"""	
		ss = Chou2006()
		ss.exptype = "fullinfo"
		ss.equations = [1]
		ar = ALRSolver(ss)
		ar.solve(l1penalty=0.0,maxiter=10000,tol=10e-6)
		assert len(ar.all_exp_art) == 1
		assert len(ar.all_exp_art[0]['eqns']) == 1
		art = ar.all_exp_art[0]['eqns'][0]			
		sol =  	{'alpha': 11.177018643003306,
		 'beta': 8.7027543289527056,
		 'g': np.array([-1.0236953]),
		 'h': np.array([ 0.62268492])}
		params = art.params[-1]
		assert np.allclose(sol['alpha'],params['alpha'],atol=1)	
		assert np.allclose(sol['beta'],params['beta'],atol=1.5)	
		assert np.allclose(sol['g'],params['g'],atol=0.5)	
		assert np.allclose(sol['h'],params['h'],atol=0.5)	

	@attr('slow')
	def test_noinfo1(self) : 
		""" Run ALR eqn 1 under no info with tol 10e-6 """
		ss = Chou2006()
		ss.exptype = "noinfo"
		ss.equations = [1]
		ar = ALRSolver(ss)
		ar.solve(l1penalty=0.0,maxiter=10000,tol=10e-6)
		assert len(ar.all_exp_art) == 1
		assert len(ar.all_exp_art[0]['eqns']) == 1
		art = ar.all_exp_art[0]['eqns'][0]		
		sol = {'alpha': 12,
			 'beta': 10,
			 'g': np.array([ 0,0,-0.8,0]),
			 'h': np.array([ 0.5,0,0,0])}
		params = art.params[-1]
		assert np.allclose(sol['alpha'],params['alpha'],atol=1.0)	
		assert np.allclose(sol['beta'],params['beta'],atol=1.0)	
		assert np.allclose(sol['g'],params['g'],atol=0.5)	
		assert np.allclose(sol['h'],params['h'],atol=0.5)	

	@attr('slow')
	def test_noinfo2(self) : 
		""" Run ALR eqn 1 under noinfo with tol 10e-7 """
		ss = Chou2006()
		ss.exptype = "noinfo"
		ss.equations = [1]
		ar = ALRSolver(ss)
		ar.solve(l1penalty=0.0,maxiter=10000,tol=10e-7)
		assert len(ar.all_exp_art) == 1
		assert len(ar.all_exp_art[0]['eqns']) == 1
		art = ar.all_exp_art[0]['eqns'][0]		
		params = art.params[-1]
		assert params is None	
				




