""" This class generates the diff eqn specified by the Chou 2006 model """

from base import SSystem

class Chou2006(object):
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

	
if __name__ == '__main__' : 
	c = Chou2006()

