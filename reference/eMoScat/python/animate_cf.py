#!/usr/bin/env python3

import pyArrays as qsc

 # First append the location of C++ library to the python path
#import sys
#sys.path.append("/opt/intel/composer/lib/intel64")
#import ctypes
#try:
#	ctypes.CDLL('libmkl_rt.so', ctypes.RTLD_GLOBAL)
#except OSError:
#	print('Intel module not found, continuing...')

#import pyArrays as em
import numpy as np
import matplotlib.pyplot as plt
from cmath import *
import json

from optparse import OptionParser

parser = OptionParser(usage="%prog [options] cf ifc tffc")
parser.add_option("-x", "--max", action='store', dest="imax", help="maximal element", metavar="int")
parser.add_option("-s", "--sampling-rate", action='store', dest="rate", help="time sampling rate [1]", metavar="rate")
parser.add_option("-e", "--elastic-channel", action='store_true', dest="elastic", help="Elastic channel switch (-1 factor in T-matrix)")

(options, args) = parser.parse_args()

smp = 1 if not options.rate else int(options.rate)
delta = 0 if not options.elastic else -1.0

if (len(args)<3):
    parser.print_usage()
    exit(-1)

src = np.loadtxt(args[0])

ifc = np.loadtxt(args[1])[499:1500:5,:]
fc = np.loadtxt(args[2])[499:1500:5,:]

print(src.shape, ifc.shape, fc.shape)

T = src[:,0]
C = src[:,1] + src[:,2]*1j
F0 = ifc[:,1] + ifc[:,2]*1j
F = fc[:,1] - fc[:,2]*1j

E = ifc[:, 0]

plt.ion()
f = plt.figure();

out = np.zeros(E.shape, np.complex)

ierg = -0.744777 #- 0.0976049
dt = T[1] - T[0]
print(dt)

Y = np.zeros([E.shape[0], T.shape[0]/smp])

val = out.copy()
er = E + ierg
cf = 1.0 / (F * F0 * 2.0 * pi)
imax = T.shape[0] if not options.imax else int(options.imax)
for i in range(0,imax):
    eit = np.exp(1j * er * T[i]) * C[i]
    val += eit * dt * cf 
    if i % smp == 0:  
        out = abs((val-delta)/ (2*pi))**2 * 4 * pi**3 / (2 * E)
        plt.clf()
        plt.plot(E,out)
        plt.draw()
        Y[:, i//smp] = out

plt.ioff()

data = { "x" : ["%.4f" % e for e in E], "y" : [["%.4f" % x for x in XX ] for XX in Y.transpose()], "frames" : imax//smp  } 
with open("animated.json", "w") as f:
    json.dump(data, f)
    f.close()

exit()
#plt.plot(E, np.real(out))
#plt.plot(E, np.imag(out))
#plt.show()

if options.save:
    path = options.save
    o = np.zeros([E.shape[0],2], np.float)
    o[:,0] = E
    o[:,1] = out
    np.savetxt(path, o)
