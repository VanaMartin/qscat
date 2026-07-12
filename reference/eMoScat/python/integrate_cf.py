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

from optparse import OptionParser

parser = OptionParser(usage="%prog [options] correlation_f init_fourier_cf tf_fourier_cf")
parser.add_option("-x", "--max", action='store', dest="imax", help="maximal element", metavar="int")
parser.add_option("-s", "--save", action='store', dest="save", help="save to file", metavar="file")
parser.add_option("-e", "--elastic-channel", action='store_true', dest="elastic", help="Elastic channel switch (-1 factor in T-matrix)")

(options, args) = parser.parse_args()

imax = 0 if not options.imax else int(options.imax)

if (len(args)<3):
    parser.print_usage()
    exit(-1)

if (imax > 0):
    src = np.loadtxt(args[0])[0:imax,:]
else:
    src = np.loadtxt(args[0])

ifc = np.loadtxt(args[1])
fc = np.loadtxt(args[2])

print(src.shape, ifc.shape, fc.shape)

T = src[:,0]
C = src[:,1] + src[:,2]*1j
F0 = ifc[:,1] + ifc[:,2]*1j
F = fc[:,1] - fc[:,2]*1j

E = ifc[:, 0]

out = np.zeros(E.shape, np.complex)

ierg = -0.0976049
#esh = ierg  + 0.125 #-0.744777
dt = T[1] - T[0]

#eit = np.exp(1j * (ierg) * T ) * C
#X = np.zeros( [T.shape[0], 3], np.float )
#X[:,0] = T
#X[:,1] = np.real(eit)
#X[:,2] = np.imag(eit)
#np.savetxt(args[0].replace("cf_", "cfie_"), X)

delta = 0 if not options.elastic else -1.0

for i in range(0, E.shape[0]):
    if abs(F[i]) > 0:
        e = E[i] + ierg;
        eit = np.exp(1j * e * T )
        val = np.inner(eit, C)
        val  = val * dt / (F[i] * F0[i]) / ( 2.0 * pi)
        out[i] = abs((val-delta)/ (2*pi))**2 * 4 * pi**3 / (2 * E[i])
    else:
        out[i] = np.nan

plt.plot(E, np.real(out))
plt.plot(E, np.imag(out))
plt.show()

if options.save:
    path = options.save
    o = np.zeros([E.shape[0],2], np.float)
    o[:,0] = E
    o[:,1] = out
    np.savetxt(path, o)
