""" A Test module to test the base module """

from nose.tools import *
from ssystem.base import SSystem, Profile, Experiment, Constraint 
import numpy.testing as npt
import numpy as npi
from ssystem.utility import *

class TestProfile(object) :
	""" Test class to test instances of Profile class"""
	def setUp(self) : 
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
	
