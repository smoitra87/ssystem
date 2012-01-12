""" Tester module for profile class """

from nose.tools import *
from ssystem.base import SSystem
from ssystem import cparser

def test_plot_profile() :
	ss = SSystem(cparser.parse('ssystem/allProblems/ss_5genes1'))	
	prof = ss.experiments[0].profile
	prof.plot_profile()

nottest(test_plot_profile)
