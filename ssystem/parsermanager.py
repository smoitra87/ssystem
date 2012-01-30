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
	def __init__(self) : 
		pass
	def setparams(self,**args) :
		pass
	def get_chou2006(self) : 
		pass
	def get_5genes1(self) : 
		fname = 'allProblems/ss_5genes1'
		ss = SSystem(cparser.parse(fname))
		modify_ss(ss) # make modifications to ss
		return ss;


def modify_ss(ss) : 
	""" Make modifications to ss according to type"""
	pass
