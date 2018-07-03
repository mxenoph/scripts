#!/usr/bin/env python

# ---------------------------------------------------------
# .npz created from count-tags.py contains 2 arrays.
# This script reads in the npz and saves the 2 arrays as .csv to be read in R

# Import modules# {{{
# Needed for division to return float
import os, sys, argparse, re
import numpy as np
# }}}

# Parsing command line arguments and creating output subdirectories# {{{
parser = argparse.ArgumentParser()
# Required arguments
parser.add_argument('-n', '--npz', nargs = '+', metavar = "file", type = str, required = True, help = 'One or more .npz files')

args = parser.parse_args()
#}}}

npz = [ f for f in args.npz if f.endswith(".npz") ]
for n in npz:
    data = np.load(n)
    for key, value in data.items():
        np.savetxt(os.path.splitext(n)[0] + '.' + key + '.csv', value)


