"""
This script is the overall manager for all the scripts that run under 
this experimental framework of s-systems. 
It makes calls to the 
parsermanager  : 
mathods : 
logger : 
author: Subhodeep Moitra (smoitra@cs.cmu.edu)
BSD License
"""

from parsermanager import ParserManager
import sys,os

class ExperimentManager(object) : 
	""" Class for managing all experiments """
	def __init__(self) : 
	
		# Create Logging directory if it does not exist
		paths = os.path.split 
		pathr = os.path.realpath
		pathd = os.path.dirname
		pathj = os.path.join
		pathe = os.path.exists
		scriptname = sys.argv[0]	
		scriptdir = pathr(pathd(scriptname))
		self.basedir = paths(scriptdir)[0]
		
		logdir = pathj(self.basedir,'logs')
		if not os.path.exists(logdir) : 
			os.mkdir(logdir)
			
		# Initialize ParserManager
		self.pman = ParserManager()
		
		# A schedule of scenarios to run
		self.schedule = [
			(self.run_scenario1,"Chou2006 experiments"),
			(self.run_scenario2,"ss_5genes1 experiments")
		]


	def run_scenario1(self) : 
		""" Runs experiments for Chou2006"""
		for ss in self.pman.get_gen_chou2006() :
			print "Running %s" %(ss.name)

	def run_scenario2(self) : 
		""" Runs experiments for ss_5genes1"""
		for ss in self.pman.get_gen_5genes1() :
			print  "Running %s" % (ss.name)
		
	def run_scenario3(self) : 
		""" Runs all experiments in all s-systems"""
		pass
	
	def set_params(self,**args) : 
		""" Sets params of the Experiment Manager Class"""
		pass

	def start(self) : 
		""" Starts execution of Exoeriment Engine"""
		# Some place for Tkinter GUI code initialization etc..
	
		
		
	
		# Execute all scenarios in schedule
		for scenario in self.schedule : 
			run_scenario,description = scenario	
			run_scenario()
			# Some code for storing results ???
			# or let scenarios handle it themselves ?

if __name__ == '__main__'  :
	eman = ExperimentManager()
	eman.start()
