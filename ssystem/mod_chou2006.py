
from utility import basedir,logdir

from other_models import Chou2006
import copy


""" Contains modifiers for Chou2006 class"""

class ModifierChou2006(object) : 
	""" Makes modifications to an s-system and returns a list of 
	modified s-systems"""
	def __init__(self) : 
		# List of mods to be applied 
		mod_list1 = [self.modNone,self.mod1,self.mod2,self.mod3,\
		self.mod4]
		mod_list2 = [self.mod6]
		mod_list3 = [self.mod4] # works
		self.mods = mod_list2
		
				

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
		""" Make algo use fullinfo
		"""
		# Set the experiment type to be under full information
		ss.exptype = "fullinfo"
		return ss

	def mod2(self,ss) : 
		""" Make algo use partial info""" 
		ss.exptype = "partialinfo"
		return ss

	def mod3(self,ss) : 
		""" Ask algo to return structure""" 
		ss.exptype = "structure"
		return ss

	def mod4(self,ss) : 
		""" Ask algo to run under full info for eqn1 only """
		ss.exptype = "fullinfo"
		ss.equations = [1]
		return ss

	def mod5(self,ss) :
		""" Ask algo to run under full info for eqn 2 only """
		ss.exptype = "fullinfo"
		ss.equations = [2]
		return ss

	def mod6(self,ss) : 
		""" Ask algo torun under no info for eqn 1 only """
		ss.exptype = "noinfo"
		ss.equations = [1]
		initsol = ss.constraint.initsol
		initsol.beta['defaultInitialValue'] = 5.0
		initsol.h['defaultInitialValue'] = 1.0
		mspace = ss.constraint.modelspace
		mspace.alpha['defaultUpperBound'] = 50.0
		mspace.beta['defaultUpperBound'] = 50.0
		mspace.g['defaultLowerBound'] = -10.0
		mspace.g['defaultUpperBound'] = 10.0
		mspace.h['defaultLowerBound'] = -10.0
		mspace.h['defaultUpperBound'] = 10.0
		return ss

	def modNone(self,ss) : 
		""" No op """
		return ss


if __name__ == '__main__' : 
	ss = Chou2006()
	modifier = ModifierChou2006()
	ss_list = modifier.gen_modifications(ss)
	

