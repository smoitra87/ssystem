""" The utility module for the ssystem package """

import numpy as np
from scipy import interpolate

dbglevel = 1


def calc_slope(profile) :
	""" 
	Accepts a profile class object and returns the slopes calculated
	by performing 
	 * spline fitting
	"""
	var = profile.var; # 2D array of biochemical profile vars
	n_sample,n_var = var.shape 
	time = profile.time; # vector of time points
	f1 = lambda(x):interpolate.splrep(time,x)
	f2 = lambda(tck):interpolate.splev(time,tck,der=1) 
	tcks = (map(f1,var.T)) # params
	#Calculate the derivatives
	derivatives = np.array(map(f2,tcks)).T
	return derivatives,tcks	

def same_dict(x,y):	
	"""
	Accepts two dicts x and y ; returns if they are the same 
	"""
	n_common = len(set(x.items()) & set(y.items()) )
	return len(x)==len(y) and n_common==len(x)
