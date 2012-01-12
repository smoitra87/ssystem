from nose.tools import *
import ssystem

def setup() :
   print "Setting up!"

def teardown() : 
   print "Tearing down!"

def test_basic() : 
   print "A test ran!"

def test_fail() :
	print "I am going to pass"
	assert True

def test_ignore() :
	pass

def Iamnotpariksha() :
	pass

def docpariksha() :
	"""
	I am a doc test
	>>> print('Hello World')
	Hello World
	"""
	pass

istest(Iamnotpariksha)
nottest(test_ignore)

