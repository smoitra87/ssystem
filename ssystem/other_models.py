""" This class generates the diff eqn specified by the Chou 2006 model """

from base import SSystem, Experiment, Profile
import numpy as np
import pylab as pl
from scipy.integrate import ode
from utility import dbglevel


class Chou2006(SSystem):
	""" This class defines the Chou 2006 model
	It is a sublass of the SSystem class
	Several experiments are generated with varying number of points, 
	noise etc.
	An SSystem model is created, which is returned through the 
	get_ssystem() method
	
	"""
	def __init__(self) : 
		ss = {
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
			'experiments' : None,
			'modelspace': self._gen_modelspace(),
			'initbound': self._gen_initbound(),
			'initsol' : self._gen_initsol()
		}
		self.ss_dict = ss
		self.ss_params = {
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
		self.y0,self.t0 = [1.4,2.7,1.2,0.4],0 
		self._calc_slope_var(0,5,50)
	

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
		for i in range(n) : 
			var = {
				'id' : i+1,
				'name' : 'x'+str(i+1),
				'property' : 'dependent'
			}	
			yield var

	def _gen_modelspace(self):
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
		initbound = {
			'alpha' : {},
			'beta' : {},
			'g' : {},
			'h' : {}
		}
		return initbound
	
	def _gen_initsol(self) :
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
	
	def _calc_slope_var(self,begin,t_end,nPoints) :
		""" Calculate the slope and X for s-system """
		dy = self._dy # set the gradient function
		y0 = self.y0 # start point
		t0 = self.t0 # start time (doesn't really matter)
		ss_params = self.ss_params # params of equations
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
		self._t = np.array(t)
		self._y = np.array(y)
		self._slope = np.array(slope)

		if(dbglevel > 2) :
			ax = pl.gca()
			ax.set_color_cycle(['r','b','g','c'])
			pl.plot(t,y,'.')
			pl.xlabel('Time')
			pl.ylabel('Concentration')
			pl.title('Plot of Conc vs time for canon Chou 2006 ssystem')
			pl.show()	
	

	def get_ssytem(self) : 
		return self.ss

	def dummy_f1(self) : 
		return self


	
if __name__ == '__main__' : 
	c = Chou2006()

