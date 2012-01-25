from nose.tools import *
from ssystem.other_models import Chou2006
from ssystem.base import Constraint
from ssystem.utility import same_dict
import numpy as np
import sys

class TestChou2006(object) :
	def setUp(self) : 
		"""Setup fixture for the test"""
		self.ss = Chou2006() # create ss instance
				
	
	def tearDown(self) :
		pass
	
	def test_keys(self) :
		""" Test Chou2006 keys"""
		ss = self.ss
		assert ss.name == 'chou2006'
		assert ss.systype == 'SSystem'
		eq_(ss.date,"01/01/12")
		eq_(ss.url,"http://www.cs.cmu.edu/~subhodee")
	
	def test_variables(self) :
		""" Test Chou2006 variables"""
		# number of vars is 5
		ss = self.ss
		var_list = [
		{'id':1,'name':'x1','property':'dependent'},
		{'id':2,'name':'x2','property':'dependent'},
		{'id':3,'name':'x3','property':'dependent'},
		{'id':4,'name':'x4','property':'dependent'}
		]
		assert len(var_list) == len(ss.variables)
		for ii in range(len(var_list)) :
			var_test = var_list[ii]
			var_ss = ss.variables[ii]
			assert len(var_ss) == len(var_test)
			eq_(var_ss['id'],var_test['id'])
			eq_(var_ss['name'],var_test['name'])
			eq_(var_ss['property'],var_test['property'])

	def test_errorfunc(self) : 
		""" Test Chou2006 errorfunc"""
		ss = self.ss
		erf = ss.errorfn
		eq_(erf['equation'],'-L+lambda*K')
		assert_almost_equal(erf['lambda'],1.0)
		eq_(erf['type'],'minusLogLikelihoodPlusLambdaK')

	def test_experiment(self) :
		""" Test some general Chou2006 exp properties"""
		ss = self.ss
		assert_greater(len(ss.experiments),0)

	def test_experiment1(self) :
		""" Test Chou2006 first experiment"""
		# Check nPoints
		ss = self.ss
		e = ss.experiments[0]
		p = e.profile
		eq_(p.n_sample,50)
		eq_(p.n_var,4)
		eq_(p.var.shape,(p.n_sample,p.n_var))
		eq_(len(p.time),p.n_sample)
		# check profile time accurate to 2 decimal places
		assert_f = lambda t : assert_almost_equal(t[0],t[1],places=2)
		map(assert_f,zip(p.time,np.arange(0,5,0.1)))
		eq_(p.time.shape,(p.n_sample,))
		# Check initial point
		assert_f = lambda t : assert_almost_equal(t[0],t[1],places=2)
		map(assert_f,zip(ss._y0,[1.4,2.7,1.2,0.4]))
		eq_(ss._t0,0)

		# check some random points accurate to 3 decimal places
		assert_f = lambda t : assert_almost_equal(t[0],t[1],places=3)
		test_vars = [
			[ 1.4       ,  2.7       ,  1.2       ,  0.4       ],#0th
			[ 0.32173054,  2.43086341,  2.6428591 ,  0.13536515], #10th
			[ 0.33980311,  1.8215696 ,  2.41620586,  0.12619968], # 20th
			[ 0.41154004,  1.95905513,  2.17694858,  0.14416145],# 30th
			[ 0.40633804,  2.03137886,  2.21211375,  0.14475924], # 40th
			[ 0.39767285,  2.01128965,  2.23608434,  0.14247661]] #49th	
		for i,j in enumerate((0,10,20,30,40,49)) : 
			map(assert_f,zip(test_vars[i],p.var[j]))

	def test_exp1_fitted1(self) : 
		"""Compares fitted vs true slopes with 1 decimal accuracy"""
		# check simulated slopes and fitted slopes are accurate to 
		# 2 decimal places
		ss = self.ss
		p = ss.experiments[0].profile
		assert_f = lambda t : assert_almost_equal(t[0],t[1],places=1)
		assert_f2 = lambda t : map(assert_f,zip(t[0],t[1]))
		map(assert_f2,zip(ss._slope,p.slopes))

		# check simulates slopes don't
	
	def test_exp1_fitted2(self) : 
		""" Compares fitted bs true slopes with 2 decimal accuracy , 
			Supposed to fail"""	
		ss = self.ss
		p = ss.experiments[0].profile
		assert_f = lambda t : assert_almost_equal(t[0],t[1],places=2)
		assert_f2 = lambda t : map(assert_f,zip(t[0],t[1]))
		
		try : 
			map(assert_f2,zip(ss._slope,p.slopes))
			assert False
		except AssertionError :
			sys.stdout.write("Supposed to fail at 2 decimal places\n")	
			pass
	


		
	def test_constraint(self) : 
		""" Test Chou2006 constraints """
		ss = self.ss
		cons = ss.constraint.modelspace
		test_cons = gen_constraint()
		test_model = test_cons.modelspace
		assert same_dict(test_model.alpha,cons.alpha)
		assert same_dict(test_model.beta,cons.beta)
		assert same_dict(test_model.g,cons.g)
		assert same_dict(test_model.h,cons.h)

		cons = ss.constraint.initsol
		test_initsol = test_cons.initsol
		assert same_dict(test_initsol.alpha,cons.alpha)
		assert same_dict(test_initsol.beta,cons.beta)
		assert same_dict(test_initsol.g,cons.g)
		assert same_dict(test_initsol.h,cons.h)
		

def gen_constraint() :
	modelspace = gen_modelspace()
	initbound = gen_initbound()
	initsol = gen_initsol()
	constraint = Constraint(modelspace,initbound,initsol)
	return constraint


def gen_modelspace():
	modelspace = {
		'alpha': {
			'defaultLowerBound' : 0.0,
			'defaultUpperBound' : 15.0
		},	
		'beta': {
			'defaultLowerBound' : 0.0,
			'defaultUpperBound' : 15.0
		},
		'g': {
			'defaultLowerBound' : -3.0,
			'defaultUpperBound' : 3.0
		},
		'h': {
			'defaultLowerBound' : -3.0,
			'defaultUpperBound' : 3.0
		}
	}
	return modelspace

def gen_initbound() :
	initbound = {
		'alpha' : {},
		'beta' : {},
		'g' : {},
		'h' : {}
	}
	return initbound

def gen_initsol() :
	initsol = {
		'alpha' : {
			'defaultInitialValue' : 1.0
		},
		'beta' : {
			'defaultInitialValue' : 1.0
		},
		'g': {
			'defaultInitiValue' : 0.0
		},
		'h' : {
			'defaultInitialValue' : 0.0
		}
	}
	return initsol

