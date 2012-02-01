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
from utility import basedir,logdir
import logging
import datetime

# This statement is needed for all future logging 
from utility import SSLogger 

class ExperimentManager(object) : 
	""" Class for managing all experiments """
	def __init__(self) : 

		# Creating logger for eman class	
		self.logger = logging.getLogger('ss.eman')
		logger = self.logger
	
		# Starting new logging sessions
		logger.info('-'*50)
		logger.info('Starting new logging session %s' % \
			datetime.datetime.now().__str__())

		# Create Logging directory if it does not exist
		logger.info('Checking if logs/ already exists') 
		if not os.path.exists(logdir) : 
			logger.info('Could not find logs/ Creating...')
			os.mkdir(logdir)
			
		# Initialize ParserManager
		logger.info('Initializing parsermanager')
		self.pman = ParserManager()
		
		# A schedule of scenarios to run
		logger.debug('Current Scenario schedule : ')
		self.schedule = [
			(self.run_scenario1,"Chou2006 experiments"),
			(self.run_scenario2,"ss_5genes1 experiments")
		]
		logger.debug([t[1] for t in self.schedule])


	def run_scenario1(self) : 
		""" Runs experiments for Chou2006"""
		logger = logging.getLogger('ss.eman.sc1')
		for ii,ss in enumerate(self.pman.get_gen_chou2006()) :
			 logger.info("Running ss: %s mod: %d" %(ss.name,ii))

	def run_scenario2(self) : 
		""" Runs experiments for ss_5genes1"""
		logger = logging.getLogger('ss.eman.sc2')
		for ii,ss in enumerate(self.pman.get_gen_5genes1()) :
			logger.info("Running ss: %s mod: %d" % (ss.name,ii))
		
	def run_scenario3(self) : 
		""" Runs all experiments in all s-systems"""
		logger = logging.getLogger('ss.eman.sc3')
		logger.info('Running all ssystems and modifications')
	
	def set_params(self,**args) : 
		""" Sets params of the Experiment Manager Class"""
		self.logger.info('Setting eman params')
		pass

	def start(self) : 
		""" Starts execution of Experiment Engine"""
		# Some place for Tkinter GUI code initialization etc..
	
		
		logger = logging.getLogger('ss.eman.start')
		
		logger.info('Starting execution of scenarios')
		# Execute all scenarios in schedule
		for scenario in self.schedule : 
			run_scenario,description = scenario	
			logger.info('Running scenario')
			run_scenario()
			# Some code for storing results ???
			# or let scenarios handle it themselves ?

if __name__ == '__main__'  :
	eman = ExperimentManager()
	eman.start()
