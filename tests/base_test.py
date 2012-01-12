""" A Test module to test the base module """

from nose.tools import *
from ssystem.base import SSystem, Profile, Experiment, Constraint 
import numpy.testing as npt
import numpy as np
from ssystem.utility import dbglevel

def gen_sample(n) :
	""" Generate n samples """
	for i in range(n) : 
		sample = {	
			'time' : i*0.1,
			'id' : i+1,
			'sdev' : [0.01]*5,
			'var' : [0.1*i,0.2*i,0.3*i,0.4*i,0.5*i]
		}
		yield sample		

def gen_profile(n) :
	""" Generate profile """
	samples = [s for s in gen_sample(n)]
	prof = Profile(samples)
	return prof

def gen_samples(n) :
	""" Generate list of samples """
	samples = [s for s in gen_sample(n)]
	return samples

def gen_experiment(n) :
	""" Generate an experiment instance """
	for i in range(n) :
		exp  = {
			'id' : i+1,
			'name' : 'exp'+str(i+1),
			'datatype' : 'perfectData',
			'samples' : gen_samples(4+i)
		}
		yield exp

def gen_vars(n)  : 
	for i in range(n) : 
		var = {
			'id' : i+1,
			'name' : 'x'+str(i+1),
			'property' : 'dependent'
		}	
		yield var

def gen_ssystem(n) :
	""" Generate an ssystem instance """
	experiments = [exp for exp in gen_experiment(n)]
	ss = {
		'name' : 'ss_dummy',
		'type' :  'SSystem',
		'date' : '01/01/12', 
		'url' : 'http://www/cs.cmu.edu/~subhodee',
		'errorfunc' : {
			'equation': '-L+lambda*K',
			'lambda': 1.0,
			'type': 'minusLogLikelihoodPlusLambdaK'			
		},
		'variables' : [var for var in gen_vars(5)],
		'experiments' : experiments,
		'constraint' : None
	}
	return ss
		

class TestProfile(object) :
	""" Test class to test instances of Profile class"""
	def setUp(self) : 
		self.samples = [s for s in gen_sample(4)]
		self.prof = Profile(self.samples)

	def tearDown(self) : 
		pass	
	
	def test_keys(self) :
		p = self.prof
		npt.assert_almost_equal(p.time,np.array([0.0,0.1,0.2,0.3]))
	
	def test_var(self) : 
		p = self.prof
		npt.assert_almost_equal(p.var,np.array(
			[[0,0,0,0,0],
			[0.1,0.2,0.3,0.4,0.5],
			[0.2,0.4,0.6,0.8,1.0],
			[0.3,0.6,0.9,1.2,1.5]
			]))

	def test_sdev(self) :
		p = self.prof
		npt.assert_almost_equal(p.sdev,np.array([[0.01]*5]*4))

	def test_slopes(self) :
		p = self.prof
		npt.assert_almost_equal(p.slopes,np.array(
			[[1,2,3,4,5]]*4),decimal=2)
		if(dbglevel > 2) : 
			print "slope shape", self.prof.slopes.shape
			self.prof.plot_profile()
	
	def set_profile(self,_prof) :
		self.prof = _prof

				

class TestExperiment(object) :
	def setUp(self) : 
		self.exp_dict = [exp for exp in gen_experiment(1)][0]
		self.exp = Experiment(self.exp_dict);

	def tearDown(self) :
		pass

	def test_keys(self) :
		exp = self.exp 
		assert_equal(exp.id,1)		
		assert_equal(exp.name,'exp1')
		assert_equal(exp.datatype,'perfectData')

	def test_samples(self) :
		exp = self.exp 
		test_prof = TestProfile() ;
		test_prof.set_profile(exp.profile)
		test_prof.test_keys()
		test_prof.test_var()
		test_prof.test_sdev()
		test_prof.test_slopes()
	
	
class TestSSystem(object) :
	def SetUp(self) : 
		pass
	
	def tearDown(self) :
		pass
	
	def test_keys(self) :
		pass


# Dummy Test Case

x = 0

def setup_fn() :
	global x
	x = 1

def teardown_fn() : 
	pass


@with_setup(setup_fn,teardown_fn)
def test_f1() : 
	eq_(x,1)
	
