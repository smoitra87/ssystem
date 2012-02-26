""" The utility module for the ssystem package """
import numpy as np
from scipy import interpolate
import logging, sys,os
import pylab as pl


dbglevel = 1
loglevel = logging.DEBUG
eps = 10e-7

def plot_pair(v1,v2,labels=['v1','v2']) :
	""" Plot a pair of vectors """
	assert len(v1) == len(v2)
	nPoints = len(v1)
	pl.plot(range(nPoints),v1,range(nPoints),v2)
	pl.legend(labels)
	pl.show()

def get_basedir() :
	""" Get the base directory of the project"""
	paths = os.path.split 
	pathr = os.path.realpath
	pathd = os.path.dirname
	pathj = os.path.join
	pathe = os.path.exists
	scriptname = sys.argv[0]	
	scriptdir = pathr(pathd(scriptname))
	basedir = paths(scriptdir)[0]
	return basedir

def get_logdir() : 
	""" Get the directory of the logs folder"""
	logdir = os.path.join(basedir,'logs')	
	return logdir

def get_homedir() : 
	""" Get the home directory of the user"""
	homedir = os.getenv('HOME')
	return homedir

basedir = get_basedir()
logdir = get_logdir()

def calc_slope(profile) :
	""" 
	Accepts a profile class object and returns the slopes calculated
	by performing 
	 * spline fitting
	"""
	logger = logging.getLogger('ss.util.spline')
	logger.debug('Estimating slopes with splines')
	var = profile.var; # 2D array of biochemical profile vars
	n_sample,n_var = var.shape 
	time = profile.time; # vector of time points
	f1 = lambda(x):interpolate.splrep(time,x)
	f2 = lambda(tck):interpolate.splev(time,tck,der=1) 
	tcks = (map(f1,var.T)) # params
	#Calculate the derivatives
	derivatives = np.array(map(f2,tcks)).T
	return derivatives,tcks	

def within_tuple(tup,x) : 
	""" Checks whether an element x is within a range specified by 
	tuple
	"""
	return x >= tup[0] and x <= tup[1]

def same_dict(x,y):	
	"""
	Accepts two dicts x and y ; returns if they are the same 
	"""
	n_common = len(set(x.items()) & set(y.items()) )
	return len(x)==len(y) and n_common==len(x)

class SSLogger(object) : 
	""" Contains Loggers, Handlers and Formatters for SSystem
Mostly contains static members that can be imported directly without reinitializing"""

	if not os.path.exists(logdir) : 
		try:
			os.mkdir(logdir)
		except OSError : # OS did not allow permission to that dir
			homedir = get_homedir() 
			logdir = os.path.join(homedir,'sslogs')
			if not os.path.exists(logdir):
				os.mkdir(logdir)
		

	_handler_args = {
		'filename' : os.path.join(logdir,'ssystem.log') ,
		'mode' : 'a'
	}
	# Handler for writing to ssystem.log
	ssh = logging.FileHandler(**_handler_args) 
	ssh.setLevel(logging.INFO)
	
	_handler_args = {
		'filename' : os.path.join(logdir,'debug.log') ,
		'mode' : 'w'
	}
	# Handler for writing to debug.log , more detailed logs
	dh = logging.FileHandler(**_handler_args) 
	dh.setLevel(logging.DEBUG)

	ch  = logging.StreamHandler() 
	ch.setLevel(logging.WARNING)

	fmt_ch_str = "%(name)s: %(levelname)s %(message)s"
	fmt_ch = logging.Formatter(fmt_ch_str)
	ch.setFormatter(fmt_ch)

	fmt_ssh_str = "%(name)s: %(levelname)s %(message)s"
	fmt_ssh = logging.Formatter(fmt_ssh_str)
	ssh.setFormatter(fmt_ssh)
	
	fmt_dh_str = "%(name)s: %(levelname)s %(message)s"
	fmt_dh = logging.Formatter(fmt_dh_str)
	dh.setFormatter(fmt_dh)

	sslogger = logging.getLogger("ss")	
	sslogger.setLevel(loglevel)
	sslogger.addHandler(ssh)
	sslogger.addHandler(dh)
	sslogger.addHandler(ch)
	#	sslogger.addHandler(dh)

	def __init__(self) : 
		pass


if __name__ == '__main__' : 
	logger = SSLogger.sslogger
	logger.info("Test")
	
	logger2 = logging.getLogger('ss.foo')
	logger2.info('Foo test')	
	logger2.debug('Foo Debug test')
	logger.warning('Serious stuff..!')

	print get_logdir()
	print get_basedir()
