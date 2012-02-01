"""
The manager for all parser and SSystem data operations

The role of this module to provide the datasets to the Experiment 
manager. Each of the s-systems need to be processed differently since
they have unique situations. This class makes sure that all the 
s-system models passed adhere to these requirements. 

It provides a generator for getting all the s-system models or
getting individual models

Author: Subhodeep Moitra(smoitra@cs.cmu.edu)
BSD License

"""

import cparser
from base import SSystem
from utility import dbglevel
import sys,os
import logging
from modifiers import *
from other_models import Chou2006

class ParserManager(object) :
	""" Manages loading, parsing and modifying S-systems"""
	def __init__(self) : 
		""" Initialize the ParserManager class"""
		self.logger = logging.getLogger('ss.pman')
	def setparams(self,**args) :
		""" Set params of ParserManager class"""
		pass
	def get_chou2006(self) : 
		""" Returns a list of Chou2006 modified models"""
		self.logger.info('Creating a canonical instance of chou2006')
		ss = Chou2006()
		ss_list = self._list_modifications(ss,ModifierChou2006())
		self.logger.debug('Returning modified list')
		return ss_list
		
	def get_5genes1(self) : 
		""" Returns a list of ss_5gene ssystems with or w/o 
		modifications"""
		fname = 'allProblems/ss_5genes1'
		self.logger.info('Creating canonical instance of ss_5genes1')
		self.logger.debug('Calling cparser on 5genes1')
		ss = SSystem(cparser.parse(fname))
		ss_list = self._list_modifications(ss,Modifier5genes1())
		self.logger.debug('Returning list of mods of 5genes1')
		return ss_list

	def get_gen_chou2006(self):
		""" Returns a generator for chou2006 with modifications"""
		self.logger.debug("Getting generator for chou2006")
		ss_list = self.get_chou2006()
		for ss in ss_list : 
			yield ss
	
	def get_gen_5genes1(self) : 
		""" Returns a generator for 5genes1 with modifications"""
		self.logger.debug("Creating generator for 5genes1")
		ss_list = self.get_5genes1() 
		for ss in ss_list : 
			yield ss
	
	def get_gen_all(self) : 
		""" Returns a generator for all experiments possible """
		self.logger.debug('Generating all experiments')
		for ss in self.get_gen_5genes1() : 
			yield ss
		
		for ss in self.get_gen_chou2006() : 
			yield ss
			
	def _list_modifications(self,ss,modObj) :
		""" Takes in a ssystem and a modifier object, makes 
modifications to the s-system and then returns a list of modified
ssystems """
		self.logger.info("Creating a list of modified ssystems")
		ss_list = modObj.gen_modifications(ss)
		return ss_list
	

if __name__ == '__main__'  :
	""" For testing Purposes"""
	pman = ParserManager()
	#ss_5genes1 = pman.get_5genes1()
	#ss_chou2006 = pman.get_chou2006()
	ss_list = [ss for ss in pman.get_gen_all() ] 
