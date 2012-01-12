#! /usr/bin/env python

"""
Base class that contains a number of important classes
author:Subhodeep Moitra(smoitra@cs.cmu.edu)

BSD License
"""

import numpy as np
import pylab as pl
from utility import calc_slope
from scipy import interpolate
import pdb, sys


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
		self.slopes,self.tcks = calc_slope(self)
		self.n_sample, self.n_var = self.var.shape
	def set_params(self,**args) :
		pass
	def plot_profile(self) : 
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
			dy = self.slopes * dx
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
	def __init__(self) : 
		self.isdefined = False;
	def set_params(self,**args) : 
		pass

