
from utility import basedir,logdir

from other_models import Chou2006
import copy


""" Contains modifiers for Chou2006 class"""

class ModifierChou2006(object) : 
	""" Makes modifications to an s-system and returns a list of 
	modified s-systems"""
	def __init__(self) : 
		# List of mods to be applied 
		self.mods = [self.modNone,self.mod1]

	def gen_modifications(self,ss)  : 
		""" Generate a list of all modifications possible to ss """
		ss_list = []
		# Cycle over the mods and make the modifications
		for mod in self.mods : 
			# Make shallow copy of ss
			ss_mod = copy.copy(ss)		
			ss_list.append(mod(ss_mod))
		return ss_list

	def mod1(self,ss) : 
		""" Throws away all but the first equation 
		Asks algorithm to be strict when dealing with this type
		"""
		# Set the experiment type to be under full information
		ss.exptype = "fullinfo"
		return ss

	def modNone(self,ss) : 
		""" No op """
		return ss


if __name__ == '__main__' : 
	ss = Chou2006()
	modifier = ModifierChou2006()
	ss_list = modifier.gen_modifications(ss)
	

