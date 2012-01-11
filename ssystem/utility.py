""" The utility module for the ssystem package """

import numpy as np
from scipy import interpolate

def calc_slope(profile) :
	""" 
	Accepts a profile class object and returns the slopes calculated
	by performing 
	 * spline fitting
	"""
	var = profile.var; # 2D array of biochemical profile vars
	n_sample,n_var = var.shape 
	time = profile.time; # vector of time points
	tcks = map(lambda(x):interpolate.splrep(time,x),var.T ) # params
	#Calculate the derivatives
	derivatives = map(lambda(tck):interpolate.splev(time,tck,der=1),tcks)
	return derivatives,tcks	

	
