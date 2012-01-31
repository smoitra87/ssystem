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



class ParserManager(object) :
	""" Manages loading, parsing and modifying S-systems"""
	def __init__(self) : 
		""" Initialize the ParserManager class"""
		self.lzmod = LazyModifier()

	def setparams(self,**args) :
		""" Set params of ParserManager class"""
		pass
	def get_chou2006(self) : 
		""" Returns a list of Chou2006 modified models"""
		pass
	def get_5genes1(self) : 
		""" Returns a list of ss_5gene ssystems with or w/o 
		modifications"""
		fname = 'allProblems/ss_5genes1'
		ss = SSystem(cparser.parse(fname))
		self.lzmod.modify(ss) # make modifications to ss
		return ss;


class LazyModifier(object) :
	""" Used for modification of s-ssystem models"""
	def __init__(self) :
		pass
	def modify(self,ss) :
		pass
	

