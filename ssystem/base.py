#! /usr/bin/env python

"""
Base class that contains a number of important classes
author:Subhodeep Moitra(smoitra@cs.cmu.edu)

BSD License
"""

import numpy as np

class SSystem(object) : 
	""" The ssystem class """
	def __init__(self,data) : 
		self.name = data['problemname']
		self.systype = data['type']
		self.date = data['date']
		self.url = data['url']
		self.variables = data['variables']
		self.errorfn = data['errorfunc']
		self.experiments = [Experiment(exp) for exp in data['experiments']]
		self.constraint = None;
	def set_params(self,**args) :
		pass

class Experiment(object) : 
	""" The experiment class """
	def __init__(self,exp) :
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
		self.time = np.array([sample['time'] for sample in samples])
		self.var = np.array([sample['var'] for sample in samples])
		self.sdev =np.array( [sample['sdev'] for sample in samples])
	def set_params(self,**args) :
		pass

class Constraint(object) : 
	"""The experiment constraint class """
	def __init__(self) : 
		self.isdefined = False;
	def set_params(self,**args) : 
		pass

