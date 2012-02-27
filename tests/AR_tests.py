"""
Test cases for Alternating Regression Module
"""

from nose.tools import *

from ssystem.AlternatingRegression import ARSolver
from ssystem.other_models import Chou2006
import numpy.testing as npt
import numpy as np
from numpy import linalg as LA

from ssystem.utility import same_dict
from ssystem.base import Profile




class TestARTracker(object) :
	""" Test cases for the AR Tracker """
	def setUp(self) : 
		pass
	
	def tearDown(self) :
		pass
	
	def test_bookkeep(self) :
		""" Test the bookkeep command """
		pass

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

	def test_regressors(self) : 
		""" Test _set_regressors code """

		# Test for fullinfo and 1st and 3rd eqn
		self.ss = Chou2006()
		self.ss.exptype = "fullinfo"
		self.ss.equations = [1,3]
		self.ar = ARSolver(self.ss)

		eq_(len(self.ar.regressors),2)		

		eq_(self.ar.regressors[0]['degrad'],[1])
		eq_(self.ar.regressors[0]['prod'],[3])
		eq_(self.ar.regressors[1]['degrad'],[3,4])
		eq_(self.ar.regressors[1]['prod'],[2])
	
		# Test for fullinfo and all eqns	
		self.ss = Chou2006()
		self.ss.exptype = "fullinfo"
		self.ar = ARSolver(self.ss)	

		eq_(self.ar.regressors[0]['degrad'],[1])
		eq_(self.ar.regressors[0]['prod'],[3])
		eq_(self.ar.regressors[1]['degrad'],[2])
		eq_(self.ar.regressors[1]['prod'],[1])
		eq_(self.ar.regressors[2]['degrad'],[3,4])
		eq_(self.ar.regressors[2]['prod'],[2])
		eq_(self.ar.regressors[3]['degrad'],[4])
		eq_(self.ar.regressors[3]['prod'],[1])

		# Test for partial info and all equations	
		self.ss = Chou2006()
		self.ss.exptype = "partialinfo"
		self.ar = ARSolver(self.ss)	

		eq_(self.ar.regressors[0]['degrad'],[1,3])
		eq_(self.ar.regressors[0]['prod'],[1,3])
		eq_(self.ar.regressors[1]['degrad'],[1,2])
		eq_(self.ar.regressors[1]['prod'],[1,2])
		eq_(self.ar.regressors[2]['degrad'],[2,3,4])
		eq_(self.ar.regressors[2]['prod'],[2,3,4])
		eq_(self.ar.regressors[3]['degrad'],[1,4])
		eq_(self.ar.regressors[3]['prod'],[1,4])

		# Test for no info and all equations
		self.ss = Chou2006()
		self.ss.exptype = "noinfo"
		self.ar = ARSolver(self.ss)	

		for ii in range(4) : 
			eq_(self.ar.regressors[ii]['degrad'],[1,2,3,4])
			eq_(self.ar.regressors[ii]['prod'],[1,2,3,4])

		# Test for no info and 1st and 3rd equations
		self.ss = Chou2006()
		self.ss.exptype = "noinfo"
		self.ss.equations = [1,3]
		self.ar = ARSolver(self.ss)	

		eq_(len(self.ar.regressors),2)

		for ii in range(2) : 
			eq_(self.ar.regressors[ii]['degrad'],[1,2,3,4])
			eq_(self.ar.regressors[ii]['prod'],[1,2,3,4])

	def test_core_init_params(self) : 
		""" Test _core_init_params"""
		self.ss = Chou2006()
		self.ss.exptype= "fullinfo"
		self.ss.equations = [1,3]
		self.ss.constraint.initsol.beta['beta_1'] = 5.0
		self.ss.constraint.initsol.h['h_1_1'] = 1.2
		self.ar = ARSolver(self.ss)
	
		(a_list,b_list,g_list,h_list) = self.ar._core_init_params()

		eq_(a_list,[1.0,1.0])
		eq_(b_list,[5.0,1.0])
		ass_f = lambda t : eq_(np.allclose(t[0],t[1]),True)
		map(ass_f,zip(h_list,[np.array([1.2]),np.array([0,0])]))
		map(ass_f,zip(g_list,[np.array([0]),np.array([0])]))


	def test_core_calc_design(self) : 
		""" Test the design matrix"""
		
		# Test under the case of fullinfo
		
		# Generate dummy data
		ss = Chou2006()
		ss.exptype = "fullinfo"
		ss.equations = [1]
		ar = ARSolver(ss)

		prof = ss.experiments[0].profile
					

		reg_p = [3]
		reg_d = [1]
		
		Lp = np.ones(prof.n_sample)
		Ld = np.ones(prof.n_sample)
		
		X_p = [np.log(prof.var[:,2])]
		X_d = [np.log(prof.var[:,0])]

		Lp = np.vstack((Lp,np.array(X_p))).T
		Ld = np.vstack((Ld,np.array(X_d))).T

		Cp = np.dot(LA.inv(np.dot(Lp.T,Lp)),Lp.T)
		Cd = np.dot(LA.inv(np.dot(Ld.T,Ld)),Ld.T)


		lp_list1,ld_list1 = [Lp],[Ld]
		cp_list1,cd_list1 = [Cp],[Cd]	
	
		# Generate code results
		lp_list2,ld_list2,cp_list2,cd_list2 = ar._core_calc_design(prof)

		eq_(len(lp_list2),1)		
		eq_(len(ld_list2),1)
		eq_(len(cp_list2),1)		
		eq_(len(cd_list2),1)
		assert np.allclose(Lp,lp_list2[0])
		assert np.allclose(Ld,ld_list2[0])
		assert np.allclose(Cp,cp_list2[0])
		assert np.allclose(Cd,cd_list2[0])

		# Test under the case of no info
		ss = Chou2006()
		ss.exptype = "noinfo"
		ss.equations = [1,2]
		ar = ARSolver(ss)

		prof = ss.experiments[0].profile
					
		lp_list1,ld_list1 = [],[]
		cp_list1,cd_list1 = [],[]	

		reg_p = [1,2,3,4]
		reg_d = [1,2,3,4]
		
		Lp = np.ones(prof.n_sample)
		Ld = np.ones(prof.n_sample)
		
		X_p = np.log(prof.var.T)
		X_d = np.log(prof.var.T)

		print "X_p.shape", X_p.shape
		print "X_d.shape", X_d.shape

		Lp = np.vstack((Lp,X_p)).T
		Ld = np.vstack((Ld,X_d)).T

		Cp = np.dot(LA.inv(np.dot(Lp.T,Lp)),Lp.T)
		Cd = np.dot(LA.inv(np.dot(Ld.T,Ld)),Ld.T)


		lp_list1,ld_list1 = [Lp,Lp],[Ld,Ld]
		cp_list1,cd_list1 = [Cp,Cp],[Cd,Cd]	
	
		# Generate code results
		lp_list2,ld_list2,cp_list2,cd_list2 = ar._core_calc_design(prof)

		eq_(len(lp_list2),2)		
		eq_(len(ld_list2),2)
		eq_(len(cp_list2),2)		
		eq_(len(cd_list2),2)
		assert np.allclose(Lp,lp_list2[0])
		assert np.allclose(Ld,ld_list2[0])
		assert np.allclose(Cp,cp_list2[0])
		assert np.allclose(Cd,cd_list2[0])
		assert np.allclose(Lp,lp_list2[1])
		assert np.allclose(Ld,ld_list2[1])
		assert np.allclose(Cp,cp_list2[1])
		assert np.allclose(Cd,cd_list2[1])

	def test_core_calc_sse(self) : 
		""" Test the sse """
		calc_sse = self.ar._core_calc_sse
		L = np.array([[1,2,3],[2,3,4]])
		b = np.array([1,2,3])
		y = np.array([12,15])
		assert_almost_equal(calc_sse(y,L,b),29)


	def test_core_phase1(self) :
		""" The phase1"""
		pass

	def test_core_phase2(self) :
		""" Test phase 2 """
		pass

	def test_core_monitor(self) :
		""" Test monitor """
		pass


	def test_fullinfo1(self) : 
		""" RUN AR (long)  eqn 1 under full info """	
		ss = Chou2006()
		ss.exptype = "fullinfo"
		ss.equations = [1]
		ar = ARSolver(ss)
		ar.solve(maxiter=10000,tol=10e-6)
		assert len(ar.all_exp_art) == 1
		assert len(ar.all_exp_art[0]['eqns']) == 1
		art = ar.all_exp_art[0]['eqns'][0]			
		sol =  	{'alpha': 11.177018643003306,
		 'beta': 8.7027543289527056,
		 'g': np.array([-1.0236953]),
		 'h': np.array([ 0.62268492])}
		assert np.allclose(sol['alpha'],art.params['alpha'],atol=1)	
		assert np.allclose(sol['beta'],art.params['beta'],atol=1.5)	
		assert np.allclose(sol['g'],art.params['g'],atol=0.5)	
		assert np.allclose(sol['h'],art.params['h'],atol=0.5)	

	def test_fullinfo2(self) : 
		""" RUN AR (long)  eqn 2 under full info """	
		ss = Chou2006()
		ss.exptype = "fullinfo"
		ss.equations = [2]
		ar = ARSolver(ss)
		ar.solve(maxiter=10000,tol=10e-6)
		assert len(ar.all_exp_art) == 1
		assert len(ar.all_exp_art[0]['eqns']) == 1
		art = ar.all_exp_art[0]['eqns'][0]			
		sol = {'alpha': 7.8801765469678156,
			 'beta': 2.8893110127454427,
			 'g': np.array([ 0.51119746]),
			 'h': np.array([ 0.76771705])}
		assert np.allclose(sol['alpha'],art.params['alpha'],atol=1)	
		assert np.allclose(sol['beta'],art.params['beta'],atol=1.5)	
		assert np.allclose(sol['g'],art.params['g'],atol=0.5)	
		assert np.allclose(sol['h'],art.params['h'],atol=0.5)	

		#assert_equal_dict(art.params,sol)
