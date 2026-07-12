#!/usr/bin/env python3

 # First append the location of C++ library to the python path
import sys
sys.path.append("/opt/intel/composer/lib/intel64")
import ctypes
try:
	ctypes.CDLL('libmkl_rt.so', ctypes.RTLD_GLOBAL)
except OSError:
	print('Intel module not found, continuing...')

import pyArrays as em
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

from optparse import OptionParser

x = [0,40]
y = [0,6]

### PARSING ARGUMENTS
parser = OptionParser(usage="%prog [options] input/dir")
parser.add_option("-t", "--test", action='store_true', dest="test", help="testing option")

(options, args) = parser.parse_args()

if len(args) < 1:
    parser.print_usage()
    print("please specify drectory")
    exit(0)

inDir = args[0]

state = input("state number -->")
fileName = "rydberg_sate_" + str(state) + ".bin"
psi = em.gVector2d()
vib = em.gVector()
if not psi.load(os.path.join(inDir, fileName)):
    print("error: wave function file could not be loaded")
    exit(0)
X,Y,Z = psi.export(x[0],x[1],200, y[0],y[1],200)
X,Y = np.meshgrid(X,Y)

R = np.zeros(Z.shape, np.complex)
vibration = input("vibration number -->")
while len(vibration)>0:
    vibName = "rydberg_sate_" + str(state) + "_vibration_" + str(vibration) + ".bin"
    if not vib.load(os.path.join(inDir, vibName)):
        print("error: vibration function file could not be loaded")
        exit(0)

    test = gVector2d(psi)
    test.element_w

    x,v = vib.export(y[0],y[1],200)

    for i in range(0,200):
        R[:,i] = Z[:,i] * v[i]

    fig = plt.figure( figsize=(7,5))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('Nuclear coordinate R ($a_0$)')
    ax.set_ylabel('Electronic coordinate r ($a_0$)')
    ax.set_zlabel('real part amplitude')

    #p1 = ax.plot_surface(X, Y, np.abs(R)**2, rstride=1, cstride=1, linewidth=0.5, alpha=0.3)
    p1 = ax.plot_wireframe(X, Y, np.abs(R)**2, linewidth=0.5, alpha=0.4)
    ax.view_init(elev=30, azim=40)
    #off = ax.get_zlim()[0]
    #cset = ax.contour(X, Y, np.real(R), offset=off, zdir='z')

    plt.tight_layout()

    plt.show()
    vibration = input("next vibration number (blank for exit) -->")
