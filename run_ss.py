#! / usr/bin/env python

""" Run the ssystem script"""

import ssystem
from ssystem import base

ss_dict = ssystem.cparser.parse("ssystem/allProblems/ss_5genes1");
ss = base.SSystem(ss_dict)

