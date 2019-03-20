#!/usr/bin/env python
#
# This is a front-end to the mac OSX utility 'opendiff'
# It allows to use opendiff to visualize SVN differences:
#
#    svn diff --diff-cmd opendiff.py
#
# Anonymous web donator

from sys import argv
from os import execlp
import re

left = ""
right = ""

argv.pop(0)

while argv:
    arg = argv.pop(0)
    if arg == "-u":
        pass
    elif arg == "-L":
        argv.pop(0)
    elif left == "":
        left = arg
    else:
        right = arg

execlp("opendiff", "opendiff", left, right)

