""" 
This module contains several classes implementing a number of Feature
Selection methods for S-systems. Currently the purpose is to service
AR like algorithms. 

The host of classes implemented here are : 
FeatSel - Basic skeleton class that needs to be derived in order to 
build a concrete Feature Selection class
LassoARFeatSel - Does Feature Selection by running Lasso in an AR like
scheme

author: Subhodeep Moitra(smoitra@cs.cmu.edu)

BSD License..!
"""

import sys, os

import numpy as np
import scipy as sci
import pylab as pl

import utility as util
from utility import overrides

class FeatSel(object) : 
	""" Init the Feat Sel class"""
	def __init__(self) : 
		self.features = None
		self.params = None

	def _preprocessor(self) : 
		""" Run preprocessor for Feature Selection, 
			This is a virtual class. Make it concrete..!
		"""
		pass

	def find_features(self) : 
		""" Find the relevant features
			This is a virtual class. Make it concrete..!
		"""
		pass

	

class LassoARFeatSel(FeatSel) : 
	""" Perform AR FeatSel by doing Lasso"""
	def __init__(self,ar,expid) : 
		self.ar = ar
		self.expid = expid
		self.paramspace = None		
		self.regressors = None
		self._preprocessor()

	@overrides(FeatSel)
	def _preprocessor(self) : 
		"""	
		Run preprocessor for Lasso AR Feature Selection 
		"""
		pass	

	@overrides(FeatSel)
	def find_features(self) : 
		""" Find the relevant features using Lasso AR """
		self.regressors = []
		for eqnid,eqn in enumerate(self.ar.equations) :
			feats = {'prod':[3],'degrad':[1,2]}
			self.regressors.append(feats)
	
		return self.regressors


if __name__ == '__main__' : 
	fsel = LassoARFeatSel()
	fsel.find_features()


	
