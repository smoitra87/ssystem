#! /usr/bin/env 

"""
Base class that contains a number of important classes
author:Subhodeep Moitra(smoitra@cs.cmu.edu)

BSD License
"""

import numpy as np
import pylab as pl
from utility import dbglevel, calc_slope
from scipy import interpolate
import pdb, sys
import logging
from utility import loglevel, SSLogger

logger = logging.getLogger('ss.base')

class SSystem(object) : 
	""" The ssystem class """
	def __init__(self,data) : 
		self.logger = logging.getLogger('ss.base.ss')
		self.logger.debug('Init SSystem class')
		self.name = data['problemname']
		self.systype = data['type']
		self.date = data['date']
		self.url = data['url']
		self.variables = data['variables']
		self.errorfn = data['errorfunc']
		self.experiments = [Experiment(exp) for exp in data['experiments']]
		modelspace = data['modelspace']
		initbound = data['initbound']
		initsol = data['initsol']
		self.constraint = Constraint(modelspace,initbound,initsol)
		self.exptype = 'noinfo'
		self.equations = range(1,len(self.variables))

	def set_params(self,**args) :
		pass

class Experiment(object) : 
	""" The experiment class """
	def __init__(self,exp) :
		self.logger = logging.getLogger('ss.base.exp')
		self.logger.debug('Initing Experiment Class')
		self.id = exp['id']
		self.name = exp['name']
		self.datatype = exp['datatype']
		self.profile = Profile(exp['samples'])
	def set_params(self,**args) : 
		pass

class Equation(object) :
	""" The equation class """
	def __init__(self) : 
		self.isdefined = False;
	def set_params(self,**args) : 
		pass

class Profile(object) : 
	""" The biochemical profile class """
	def __init__(self,samples) : 
		self.logger = logging.getLogger('ss.base.pro')
		self.logger.debug('Initing Profile Class')
		self.time = np.array([sample['time'] for sample in samples])
		self.var = np.array([sample['var'] for sample in samples])
		self.sdev =np.array( [sample['sdev'] for sample in samples])
		self.logger.debug('Calculating slopes')
		self.slopes,self.tcks = calc_slope(self)
		self.n_sample, self.n_var = self.var.shape
	def set_params(self,**args) :
		pass
	def plot_profile(self) :
		self.logger.debug('Plotting Profile curves') 
		try :
			xnew = np.linspace(self.time[0],self.time[-1],num=50)
			#Function handle for generating plot interpolate points
			f1 = lambda(tck) : interpolate.splev(xnew,tck)
			ynew = np.array(map(f1,self.tcks)).T
			ax = pl.gca()
			color_list = ['r','g','b','c','y']
			ax.set_color_cycle(color_list)
			pl.plot(self.time,self.var,'x',xnew,ynew);
			# plot arrows 
			arrow_scale = 40.0
			dx = (self.time[-1]-self.time[0])/arrow_scale;
			dy = self.slopes.T * dx  # transpose operation here
			if(dbglevel > 2) : 
				print "slopes shape", self.slopes.shape
			for ii in range(self.n_sample) : 
				t = self.time[ii]
				X = self.var[ii,:]
				DY = dy[:,ii]
				for jj in range(self.n_var) :
					pl.arrow(t,X[jj],dx,DY[jj],linewidth=1)
			pl.title('Plot of spline fitted biochem profile vs time')
			pl.xlabel('time')
			pl.ylabel('profile')
			pl.axis('tight')
			pl.show()
		except AttributeError :
			sys.stderr.writelines(
			"Define Profile before trying to plot..!")
			raise

class Constraint(object) : 
	"""The experiment constraint class """
	def __init__(self,modelspace,initbound,initsol) : 
		self.logger = logging.getLogger('self.base.con')
		self.logger.debug('Initing Constraint Class')
		self.modelspace = ModelSpace(modelspace)
		self.initbound = InitBound(initbound)
		self.initsol = InitSol(initsol)
	def set_params(self,**args) : 
		pass


class ModelSpace(object) :
	""" The Model space class """
	def __init__(self,modelspace):
		self.logger = logging.getLogger('self.base.mod')
		self.logger.debug('Initing Modelspace class')
		self.alpha = modelspace['alpha'];
		self.beta = modelspace['beta'];
		self.g = modelspace['g']
		self.h = modelspace['h']
	def set_params(self,**args) : 
		pass

class InitBound(object) :
	""" The Initial Bounds on variables class """
	def __init__(self,initbound) :
		self.logger = logging.getLogger('ss.base.ib')
		self.logger.debug('Initing InitBound')
		self.alpha = initbound['alpha']
		self.beta = initbound['beta']
		self.g = initbound['g']
		self.h = initbound['h']
	def set_params(self,**args) : 
		pass

class InitSol(object) :
	""" The Initial Solution class """
	def __init__(self,initsol) : 
		self.logger = logging.getLogger('ss.base.is')
		self.logger.debug('Init InitSol')
		self.alpha = initsol['alpha']
		self.beta = initsol['beta']
		self.g = initsol['g']
		self.h = initsol['h']
	def set_params(self,**args) : 
		pass


class Result(object) :
	""" Class for storing results from a single experiment"""
	def __init__(self,ss,res_method) : 
		"""
			ss - ssystem class for which the experiment was run
			res_method - result returned from the method 
		"""
		self.logger = logging.getLogger('ss.base.res')
		self.logger.debug('Initing Result Class')

	def gen_str(self) : 
		""" Generates a string representation of the results"""
		pass

	def gen_xml(self) :
		""" Generates an xml representation of the results"""
		pass

class ResultsScenario(object) : 
	""" Class for storing results from a scenario """
	def __init__(self,list_res) : 
		"""
		list_res - Contains a list of Result objects generated from various runs	 
		"""
		self.logger = logging.getLogger('ss.base.rscn')
		self.logger.debug('Initing ResultsScenario Class')

	def dump_results(self) : 
		""" Dumps all results to timestamped folder"""
		pass
	
	def _dump_pickle(self) : 
		""" Pickes all results together and dumps the results"""
		pass
	
	def _dump_xml(self) : 
		""" Dumps an xml version of results"""
		pass

	def _dump_str(self) : 
		""" Dumps a File containing results in readable format"""
		pass

	def _dump_html(self) :
		""" Dumps a html version of the results"""
		pass
