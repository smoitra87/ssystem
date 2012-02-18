""" This class generates the diff eqn specified by the Chou 2006 model """

from base import SSystem, Experiment, Profile
import numpy as np
import pylab as pl
from scipy.integrate import ode
from utility import dbglevel, loglevel
import logging

logger = logging.getLogger('ss.oth')


class Chou2006(SSystem):
	""" This class defines the Chou 2006 model
	It is a sublass of the SSystem class
	Several experiments are generated with varying number of points, 
	noise etc.
	An SSystem model is created, which is returned through the 
	get_ssystem() method
	
	"""
	def __init__(self) :
		""" Initialize the Chou 2006 class """ 
		self.logger = logging.getLogger('ss.oth.2k6') 
		self.logger.debug('Initing Chou2006 canon class')
		# True params of the S-system
		self._ss_params = {
		'alpha' : [12,8,3,2],
		'beta' : [10,3,5,6],
		'g' : [
		    [0,0,-0.8,0],
		    [0.5,0,0,0],
		    [0,0.75,0,0],
		    [0.5,0,0,0]
		    ],
		'h' : [
		    [0.5,0,0,0],
		    [0,0.75,0,0],
		    [0,0,0.5,0.2],
		    [0,0,0,0.8]
		]
		}		
		self.logger.debug('Simulating Chou2006 with Integrator')
		# Initial points of the s-system
		self._y0,self._t0 = [1.4,2.7,1.2,0.4],0 
		self._t,self._y,self._slope = self._calc_slope_var(0,5,50)
		# define the dict for the s-system
		_ss_dict = {
			'problemname' : 'chou2006',
			'type' :  'SSystem',
			'date' : '01/01/12', 
			'url' : 'http://www.cs.cmu.edu/~subhodee',
			'errorfunc' : {
				'equation': '-L+lambda*K',
				'lambda': 1.0,
				'type': 'minusLogLikelihoodPlusLambdaK'			
			},
			'variables' : [var for var in self._gen_variables(4)],
			'experiments' : self._gen_experiments(),
			'modelspace': self._gen_modelspace(),
			'initbound': self._gen_initbound(),
			'initsol' : self._gen_initsol()
		}
		self._ss_dict = _ss_dict
		# Pass ss_dict to parent class for generating s-system class
		self.logger.debug('Calling Parent SSystem constructor')
		super(Chou2006,self).__init__(_ss_dict)

	def _dy(self,t,y,ss) :
	    slope = []
	    for eqn in zip(ss['alpha'],ss['beta'],ss['g'],ss['h']) : 
	        a,b,g,h = eqn
	        expf = lambda(x) : x[0] ** x[1]
	        mulf = lambda x,y  : x*y
	        prod = a*reduce(mulf,map(expf,zip(y,g)))
	        degrad = b*reduce(mulf,map(expf,zip(y,h))) 
	        slope.append(prod-degrad)
	    return slope
	
	def _gen_variables(self,n)  :
		self.logger.debug('Create Generator for variables')
		for i in range(n) : 
			var = {
				'id' : i+1,
				'name' : 'x'+str(i+1),
				'property' : 'dependent'
			}	
			yield var

	def _gen_modelspace(self):
		self.logger.debug('Create modelspace')
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
	def _gen_initbound(self) :
		self.logger.debug('Generate Init Bounds')
		initbound = {
			'alpha' : {},
			'beta' : {},
			'g' : {},
			'h' : {}
		}
		return initbound
	
	def _gen_initsol(self) :		
		self.logger.debug('Generate InitSol')
		initsol = {
			'alpha' : {
				'defaultInitialValue' : 1.0
			},
			'beta' : {
				'defaultInitialValue' : 1.0
			},
			'g': {
				'defaultInitialValue' : 0.0
			},
			'h' : {
				'defaultInitialValue' : 0.0
			}
		}
		return initsol
	
	def _calc_slope_var(self,begin,t_end,nPoints) :
		""" Calculate the slope and X for s-system """
		dy = self._dy # set the gradient function
		y0 = self._y0 # start point
		t0 = self._t0 # start time (doesn't really matter)
		ss_params = self._ss_params # params of equations
		r = ode(dy).set_integrator('vode',method='bdf',order=15)
		r.set_initial_value(y0,t0).set_f_params(ss_params)
		
		dt = (t_end+0.) / nPoints
		t,y = [t0],[y0]
		slope = [dy(0,y0,ss_params)]
		
		# The 2*dt is to make sure that only nPoints are generated
		while r.successful() and r.t < t_end-2*dt : 
			r.integrate(r.t+dt)
			#print r.t,r.y
			t.append(r.t)
			y.append(r.y)
			slope.append(dy(r.t,r.y,ss_params))
		
		# Define the coeffs		
		_t = np.array(t)
		_y = np.array(y)
		_slope = np.array(slope)

		if(dbglevel > 2) :
			ax = pl.gca()
			ax.set_color_cycle(['r','b','g','c'])
			pl.plot(t,y,'.')
			pl.xlabel('Time')
			pl.ylabel('Concentration')
			pl.title('Plot of Conc vs time for canon Chou 2006 ssystem')
			pl.show()	

		# return time points, vars and true slope
		return _t,_y,_slope 
	
	
	def _gen_experiment(self) :
		""" Generate a single experiment fom diffeqn integration""" 
		self.logger.debug('Generate an Experiment')
		exp  = {
			'id' : 1,
			'name' : 'exp'+str(1),
			'datatype' : 'perfectData',
			'samples' : self._gen_samples()
		}
		return exp
	
	def _gen_experiments(self) : 
		""" Packages all the experiments for this class in a list"""
		self.logger.debug('Package all experiments into list')
		exp1 = self._gen_experiment()
		exps = [exp1]
		return exps

	def _gen_samples(self) :
		self.logger.debug('Generate samples')
		nSamples,nVar = self._y.shape
		samples = []
		for i,var in enumerate(self._y) : 
			sample  = {
			'time' : self._t[i],
			'id' :  str(i+1),
			'sdev' :  [0.01]*nVar,
			'var' : var
			}
			samples.append(sample)
		return samples
	
if __name__ == '__main__' : 
	c = Chou2006()

