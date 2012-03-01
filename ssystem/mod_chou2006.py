
from utility import basedir,logdir

from other_models import Chou2006
import copy


""" Contains modifiers for Chou2006 class"""

class ModifierChou2006(object) : 
	""" Makes modifications to an s-system and returns a list of 
	modified s-systems"""
	def __init__(self) : 
		# List of mods to be applied 
		#mod_list1 = [self.modNone,self.mod1,self.mod2,self.mod3,\
		#self.mod4]

		# No info cases
		mod_list2 = [self.mod7]

		# Partial info cases
		mod_list3 = [self.mod8,self.mod2]

		# Full info cases
		mod_list3 = [self.mod4_1] # works

		# all the full info cases seeded at close to some values
		mod_list4 = [self.mod4_1,self.mod4_2,self.mod4_3,self.mod4_4]
		#mod_list4 = [self.mod4_3, self.mod4_1,self.mod4_2,self.mod4_4]

		# No information case
		#mod_list5 = [self.mod5_1,self.mod5_2,self.mod5_3,self.mod5_4]
		mod_list5 = [self.mod5_1]

		self.mods = mod_list5
		
				

	def gen_modifications(self,ss)  : 
		""" Generate a list of all modifications possible to ss """
		ss_list = []
		# Cycle over the mods and make the modifications
		for mod in self.mods : 
			# Make shallow copy of ss
			ss_mod = copy.copy(ss)
			ss_mod.constraint = copy.deepcopy(ss.constraint)
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

	def mod4_1(self,ss) : 
		""" Ask algo to run under full info for eqn1 only """
		ss.exptype = "fullinfo"
		ss.equations = [1]
		return ss

	def mod4_2(self,ss) :
		""" Ask algo to run under full info for eqn 2 only """
		ss.exptype = "fullinfo"
		ss.equations = [2]
		return ss

	def mod4_3(self,ss) :
		""" Ask algo to run under full info for eqn 3 only """
		ss.exptype = "fullinfo"
		ss.equations = [3]
		initsol = ss.constraint.initsol
		initsol.beta['beta_3'] = 5.0
		initsol.h['h_3_1'] = 0.1
		initsol.h['h_3_2'] = 0.1
		initsol.h['h_3_3'] = 0.6
		initsol.h['h_3_4'] = 0.23
		return ss

	def mod4_4(self,ss) :
		""" Ask algo to run under full info for eqn 4 only """
		ss.exptype = "fullinfo"
		ss.equations = [4]
		initsol = ss.constraint.initsol
		initsol.beta['beta_4'] = 6.1
		initsol.h['h_4_1'] = 0.01
		initsol.h['h_4_2'] = 0.1
		initsol.h['h_4_3'] = 0.1
		initsol.h['h_4_4'] = 0.79

		return ss


	def mod5_1(self,ss) : 
		""" Ask algo torun under no info for eqn 1 only """
		ss.exptype = "noinfo"
		ss.equations = [1]
		return ss

	def mod5_2(self,ss) : 
		""" Ask algo torun under no info for eqn 1 only """
		ss.exptype = "noinfo"
		ss.equations = [2]
		return ss

	def mod5_3(self,ss) : 
		""" Ask algo torun under no info for eqn 1 only """
		ss.exptype = "noinfo"
		ss.equations = [3]
		return ss

	def mod5_4(self,ss) : 
		""" Ask algo torun under no info for eqn 1 only """
		ss.exptype = "noinfo"
		ss.equations = [4]
		return ss


	def mod7(self,ss) : 
		""" Ask algo to run no info for all equations """
		ss.exptype = "noinfo"
		return ss

	def mod8(self,ss) :
		""" Ask algo to run under full info for eqn 2 only """
		ss.exptype = "partialinfo"
		ss.equations = [1]
		return ss


	def modNone(self,ss) : 
		""" No op """
		return ss


if __name__ == '__main__' : 
	ss = Chou2006()
	modifier = ModifierChou2006()
	ss_list = modifier.gen_modifications(ss)
	

