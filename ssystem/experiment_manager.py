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
	def run_scenario1(self) : 
		""" Runs experiments for Chou2006"""
		pass
	def run_scenario2(self) : 
		""" Runs experiments for ss_5genes1"""
		ss = self.pman.get_5genes1()
		
	def run_all(self) : 
		""" Runs all experiments in all s-systems"""
		pass
	
	def set_params(self,**args) : 
		""" Sets params of the Experiment Manager Class"""
		pass

if __name__ == '__main__'  :
	e = ExperimentManager()
	e.run_scenario2()
