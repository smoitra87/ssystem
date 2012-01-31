""" Contains modifiers for Chou2006 class"""

class Chou2006Modifier(object) : 
	""" Makes modifications to an s-system and returns a list of 
	modified s-systems"""
	def __init__(self) : 
		pass

	def gen_modifications(self,ss)  : 
		""" Generate a list of all modifications possible to ss """
		return [ss]
	
