#!/usr/bin/env python3

# First append the location of C++ library to the python path
import sys, os
sys.path.append("/opt/intel/composer/lib/intel64")
import ctypes
try:
	ctypes.CDLL('libmkl_rt.so', ctypes.RTLD_GLOBAL)
except OSError:
	print('Intel module not found, continuing...')

import pyArrays as pa
import numpy as np
from optparse import OptionParser
import json

parser = OptionParser(usage="%prog [options] file.bin")
parser.add_option("-m", "--mode", action='store', dest="mode", help="Complex numbers interpretation [abs^2], real, imag, abs", metavar="modestring")
parser.add_option("-t", "--title", action='store', dest="title", help="sets figure title", metavar="string")

(options, args) = parser.parse_args()

x = [0, 20]
y = [0, 5]

sampling = 50

title = "Made by qscat for Plotly.js" if not options.title else options.title
mode = "abs^2" if not options.mode else options.mode

### PARSING ARGUMENTS
if len(args) < 2:
    parser.print_usage()
    exit(-1)

fname = args[0]
jname = args[1]

psi = pa.gVector2d()
if not psi.load(fname):
    print("Cannot read file", fname)
    exit(-1)

X,Y,Z = psi.export(x[0],x[1],sampling, y[0],y[1],sampling)

out = {}

if mode == "abs":
    Z = abs(Z)
elif mode=="real":
    Z = np.real(Z)
elif mode=="imag":
    Z = np.imag(Z)
else:
    Z = abs(Z)**2

out["data"] = [ {
                    "x" : ["%.4e" % x for x in X],
                    "y" : ["%.4e" % y for y in Y],
                    "z" : [ ["%.4e" % z for z in ZZ] for ZZ in Z],
                    "type" : "surface",
              } ]
out["layout"] = {"width" : 450, "height" : 500, "margin" : 0 }


with open(jname, 'w') as f:
    json.dump(out,f)
    f.close()
