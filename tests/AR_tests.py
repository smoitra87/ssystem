"""
Test cases for Alternating Regression Module
"""

from nose.tools import *

from ssystem.AlternatingRegression import ARSolver
from ssystem.other_models import Chou2006
import numpy.testing as npt
import numpy as np

from ssystem.utility import same_dict

class TestARSolverChou2006(object) :
	""" Test cases for the AR solver """
	def setUp(self) : 
		self.ss = Chou2006()
		self.ar = ARSolver(self.ss)
	
	def tearDown(self) :
		pass
	
	def test_parse_initsol(self) :
		""" Test the _parse_initsol command """
		self.ss = Chou2006()
		self.ar = ARSolver(self.ss)
		initsol = self.ss.constraint.initsol

		# Test the alphas
		initsol.alpha['defaultInitialValue'] = 0.5
		self.ar._parse_initsol()
		assert_f =  lambda t : assert_almost_equal(t,0.5,places=2)
		map(assert_f,self.ar.initsol['alpha'])


		initsol.alpha['alpha_1'] = 5.0
		initsol.alpha['alpha_3'] = 20.0
		self.ar._parse_initsol()
		assert_almost_equal(self.ar.initsol['alpha'][0],5.0,places=2)
		assert_almost_equal(self.ar.initsol['alpha'][1],0.5,places=2)
		assert_almost_equal(self.ar.initsol['alpha'][2],20.0,places=2)
		assert_almost_equal(self.ar.initsol['alpha'][3],0.5,places=2)

		# Test the betas
		initsol.beta['defaultInitialValue'] = 0.5
		self.ar._parse_initsol()
		assert_f =  lambda t : assert_almost_equal(t,0.5,places=2)
		map(assert_f,self.ar.initsol['beta'])


		initsol.beta['beta_1'] = 5.0
		initsol.beta['beta_3'] = 20.0
		self.ar._parse_initsol()
		assert_almost_equal(self.ar.initsol['beta'][0],5.0,places=2)
		assert_almost_equal(self.ar.initsol['beta'][1],0.5,places=2)
		assert_almost_equal(self.ar.initsol['beta'][2],20.0,places=2)
		assert_almost_equal(self.ar.initsol['beta'][3],0.5,places=2)
	
		# Test g
		initsol.g['defaultInitialValue'] = 0.5
		self.ar._parse_initsol()
		assert_f =  lambda t : assert_almost_equal(t,0.5,places=2)
		map(assert_f,self.ar.initsol['g'].ravel())

		initsol.g['g_1_1'] = 0.45
		initsol.g['g_3_2'] = 0.33
		self.ar._parse_initsol()
		assert_almost_equal(self.ar.initsol['g'][0][0],0.45,places=2)
		assert_almost_equal(self.ar.initsol['g'][2][1],0.33,places=2)
		eq_(self.ar.initsol['g'][1][0],0.5)
		eq_(self.ar.initsol['g'][3][1],0.5)

	def test_parse_modelspace(self) :
		""" Test the _parse_modelspace command """
		self.ss = Chou2006()
		self.ar = ARSolver(self.ss)
		modelspace = self.ss.constraint.modelspace


		# Test the alphas
		modelspace.alpha['defaultLowerBound'] = 0.5
		modelspace.alpha['defaultUpperBound'] = 30.0

		self.ar._parse_modelspace()
		assert_f =  lambda t : eq_(t,(0.5,30.0))
		map(assert_f,self.ar.modelspace['alpha'])

		modelspace.alpha['alpha_1'] = 'nonZero'
		modelspace.alpha['alpha_3'] = (1.0,20.0)
		self.ar._parse_modelspace()
		eq_(self.ar.modelspace['alpha'][0],'nonZero')
		eq_(self.ar.modelspace['alpha'][1],(0.5,30.0))
		eq_(self.ar.modelspace['alpha'][2],(1.0,20.0))
		eq_(self.ar.modelspace['alpha'][3],(0.5,30.0))

		# Test the betas
		modelspace.beta['defaultLowerBound'] = 0.5
		modelspace.beta['defaultUpperBound'] = 30.0

		self.ar._parse_modelspace()
		assert_f =  lambda t : eq_(t,(0.5,30.0))
		map(assert_f,self.ar.modelspace['beta'])

		modelspace.beta['beta_1'] = 'nonZero'
		modelspace.beta['beta_3'] = 'nonZero'
		self.ar._parse_modelspace()
		eq_(self.ar.modelspace['beta'][0],'nonZero')
		eq_(self.ar.modelspace['beta'][1],(0.5,30.0))
		eq_(self.ar.modelspace['beta'][2],'nonZero')
		eq_(self.ar.modelspace['beta'][3],(0.5,30.0))
			
		# Test g
		modelspace.g['defaultLowerBound'] = 0.5
		modelspace.g['defaultUpperBound'] = 30.0
		self.ar._parse_modelspace()
		assert_f =  lambda t : eq_(t,(0.5,30.0))
		assert_f2 = lambda t : map(assert_f,t)
		map(assert_f2,self.ar.modelspace['g'])

		modelspace.g['g_1_1'] = 'nonZero'
		modelspace.g['g_3_2'] = (0.33,20.0)
		self.ar._parse_modelspace()
		eq_(self.ar.modelspace['g'][0][0],'nonZero')
		eq_(self.ar.modelspace['g'][1][0],(0.5,30.0))
		eq_(self.ar.modelspace['g'][2][1],(0.33,20.0))
		eq_(self.ar.modelspace['g'][3][1],(0.5,30.0))
	
