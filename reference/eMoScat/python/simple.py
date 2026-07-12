#!/usr/bin/env python

 # First append the location of C++ library to the python path
import sys 
sys.path.append("/opt/intel/composer/lib/intel64")	
import ctypes
try:
	ctypes.CDLL('libmkl_rt.so', ctypes.RTLD_GLOBAL)
except OSError:
	print 'Intel module not found, continuing...'

import pyArrays as em
import numpy as np
from math import *
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import special
import time

inDir='input/tests/'
chbcf=np.zeros([1],np.complex64)
setch=False
eMax=300.0


def setCh(v, H, o, dt):
    global setch
    global chbcf
    global eMax
    if setch:
        return
    print "Setting Chebyshev approximation coefficients"
    nb = v.size()
    chbcf=np.zeros([o],np.complex64)
    # determine the spectral radius
    #W, V = H.eigenSystem().export()
    #R = -1j*dt*(W[nb-1] - W[0])/2.0
    #G = -1j*dt*W[0]
    R = -1j * dt * (eMax) / 2.0 
    G = 0 #-1j * dt * 0.154
    r = abs(R)
    for n in range(0,o):
        bess = em.bessel_jn(complex(r),n)
        #bess = special.jn(n,r)
        q = 1.0 if n==0 else 0.0
        chbcf[n] = (-1j)**n * cmath.exp(R+G) * (2.0 - q) * bess
    setch=True
    print "Chebyshev coefficients set"

# Evolution operator
def ChebEvolution(v, H, o, dt):
    global chbcf
    setCh(v,H,o,dt)
    a1=em.gVector(v)
    a2=em.gVector(v)
    out=em.gVector(v)

    if o>1000 or dt >1.0:
        print "Error! The parameters of chebyshev step exceeded possible boundaries"
        exit(1)
    
    for k in range(0,o):
        a = complex(chbcf[k])
        if k==0:
            out *= a 
        elif k==1:
            H.gemv(complex(1),v,complex(0),a1)
            out.axpy(a,a1)
            v.copy(a1)
        else: 
            H.gemv(complex(2),a1,complex(-1),a2)
            v.swap(a2)
            a1.copy(v)
            out.axpy(a,a1)
    v.swap(out)
    #return out

mgp = em.multiGridParameters(inDir+'grids.txt')
gp = mgp[0]

gX = em.grid(gp)
nbX = gX.nb()

v = em.gVector(gX)
p = em.gVector(gX)

sig = 4.0

for i in range(0,nbX):
    val = 1.0/sqrt(sqrt(sig**2*pi)) * cmath.exp(-(gX.xr(i)-40.0)**2/(2.0*sig**2) - 1.0j*gX.xr(i))
    v.setVal(val, i)
    if gX.xr(i) > 10 and gX.xr(i) < 20:
        p.setVal(0.5, i)
    else:
        p.setVal(0.0, i)

print v*v

X,F=v.export(0.0, 60.0, 501)
plt.ion()
plt.plot(X,F.real)
plt.plot(X,F.imag)
plt.draw()

Opf = em.fOperator(gX, nbX, 0)
Opf.addKineticTerm(1.0)
W, V = Opf.eigenSystem().export()
print W[0], W[nbX-1]
eMax = W[nbX-1].real


Op = em.rcOperator(gX, nbX, 0)
Op.setKineticTerm(1.0)
Op+=p
Op*=1.0/eMax

for i in range(0,1000):
    for j in range(0,10):
        ChebEvolution(v,Op,50,0.05)
    print v*v

    X,F=v.export(0.0, 60.0, 401)
    plt.gcf().clear()
    plt.plot(X,abs(F)**2)
    plt.draw()

plt.ioff()
plt.show()
