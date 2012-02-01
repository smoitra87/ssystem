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

class ExperimentManager(object) : 
	""" Class for managing all experiments """
	def __init__(self) : 
		self.pman = ParserManager()
		self.schedule = self.scenario_schedule()
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
	
	def scenario_schedule(self) : 
		"""Gives a schedule of scenarios to run """
		sched = [
			self.run_scenario1,
			self.run_scenario2
		]
		
		return sched

	def set_params(self,**args) : 
		""" Sets params of the Experiment Manager Class"""
		pass

	def start(self) : 
		""" Starts execution of Exoeriment Engine"""

		for run_scenario in self.schedule : 
			run_scenario()

if __name__ == '__main__'  :
	e = ExperimentManager()
	e.start()
