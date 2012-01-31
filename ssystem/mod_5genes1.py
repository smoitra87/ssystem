""" Contains modifiers for 5genes1 class"""

class Modifier5genes1(object) : 
	""" Makes modifications to an ss_5genes1 and returns a list of 
	modified s-systems"""
	def __init__(self) : 
		pass

	def gen_modifications(self,ss)  : 
		""" Generate a list of all modifications possible to ss """
		return [ss]
	
