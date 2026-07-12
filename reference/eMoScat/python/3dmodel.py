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

emax = 0
emin = 0

def ErgRange(gXY, A):
    global emin
    global emax

    u = em.gVector2d(gXY)
    for i in range(0, gXY.nb()):
        u[i] = complex(1.0)
    u *= complex(1.0)/ cmath.sqrt(u*u)
    v = em.gVector2d(u)
    w = em.gVector2d(u)

    alpha = complex(0)

    loop = 20000

    for i in range(0, loop):
        # normalized vector and eigenvalue (approx)
        A.gemv(complex(1), u, complex(0), v)    # Ax
        alpha = v*v                             # lambda
        v *= complex(1.0)/cmath.sqrt(alpha)     # v*v = 1.0
        u.swap(v)

    print cmath.sqrt(alpha)
    emax = cmath.sqrt(alpha)

    p = em.zVector(gXY.nb())
    
    for i in range(0, gXY.nb()):
        p[i] = -emax
  
    A += p

    for i in range(0, loop):
        # normalized vector and eigenvalue (approx)
        A.gemv(complex(1), w, complex(0), v)    # Ax
        alpha = v*v                             # lambda
        v *= complex(1.0)/cmath.sqrt(alpha)     # v*v = 1.0
        w.swap(v)

    print emax - cmath.sqrt(alpha)
    emin = emax - cmath.sqrt(alpha)


    return w

def setCh(v, H, o, dt):
    global setch
    global chbcf
    global emax
    global emin
    if setch:
        return
    print "Setting Chebyshev approximation coefficients"
    #nb = v.size()
    chbcf=np.zeros([o],np.complex64)
    # determine the spectral radius
    R = -1j * dt * (emax - emin) / 2.0 
    G = -1j * dt * emin
    r = abs(R)
    for n in range(0,o):
        #bess = em.bessel_jn(complex(r),n)
        bess = special.jn(n,r)
        q = 1.0 if n==0 else 0.0
        chbcf[n] = (-1j)**n * cmath.exp(R+G) * (2.0 - q) * bess
    setch=True

# Evolution operator
def ChebEvolution(v, H, o, dt):
    global chbcf
    setCh(v,H,o,dt)
    a1=em.gVector2d(v)
    a2=em.gVector2d(v)
    out=em.gVector2d(v)

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

def ChebEvolutionGPU(v, H, o, dt):
    global chbcf
    setCh(v,H,o,dt)
    a1 = em.CudaZVector(v)
    a2 = em.CudaZVector(v)
    out = em.CudaZVector(v)

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


mgp = em.multiGridParameters(inDir+'grids.txt')
gp = mgp[0]

gX = em.grid(gp)
gY = em.grid(mgp[1])
nbx = gX.nb()
nby = gY.nb()
nb = nbx*nby

v = em.gVector(gX)
p = em.gVector(gX)

sig = 6.0

gXY=em.grid2d(gX,gY)
gV = em.gVector2d(gXY)

# test input
R = gY.xpos() - gY.xneg() 
V = np.zeros([nbx,nby], np.complex128)
for i in range(0,nbx):
    for j in range(0,nby):
        #V[i,j] = 1.0/sqrt(sqrt(sig**2*pi)) * cmath.exp(-(gX.xr(i)-40.0)**2/(2.0*sig**2) + 0.5j*gX.xr(i)) * ( cmath.sin(gY.xr(j)*pi/R)*pi/R )
        V[i,j] = 1.0/sqrt(sqrt(sig**2*pi))**2 * cmath.exp(-(gX.xr(i)-40.0)**2/(2.0*sig**2) -(gY.xr(j)-20.0)**2/(2.0*sig**2) + 0.5j*gX.xr(i)) 
        #gV[i*nbx + j] = 1.0/sqrt(sqrt(sig**2*pi)) * cmath.exp(-(gX.xr(j)-40.0)**2/(2.0*sig**2) - 0.0j*gX.xr(j)) * ( cmath.sin(gY.xr(i)*pi)*pi )
        #V[i,j] = -2.0 * exp( - (gX.xr(i) - 10)**2/16 - (gY.xr(j) - 10)**2/16 )


gV.set(V)
gV *= 1.0 / cmath.sqrt(gV*gV)
X, Y, F = gV.export(0.0, 90.0, 50, 0.0, 40.0, 50)

x, y = np.meshgrid(X, Y)

print x.shape, y.shape, F.shape

#print "Plotting:"
#f = plt.figure()
#ax = f.add_subplot(111, projection='3d')
#ax.plot_wireframe(x, y, F)
#plt.show()

Op = em.operator2d(gXY)
Op.setKineticTerm(1.0,1.0)

Rp = em.operator2d(Op)
w = ErgRange(gXY,Rp)

R = (emax - emin).real
f = em.zVector(gXY.nb())
for i in range(0, gXY.nb()):
    f[i] = -(emax+emin).real

Op*=2.0
Op+=f
Op*= 1.0/R

gpV = em.buildGpuRep(gV.body())
gpOp = em.buildGpuRep(Op.body())
gW = em.gVector2d(gV)

order = 50
dt = 0.01

plt.ion()
t0 = time.time()
for i in range(0,10):
    for j in range(0,100):
        ChebEvolution(gV, Op, order, dt)
        #ChebEvolutionGPU(gpV, gpOp, order, dt)
   
    em.cudaPull(gW, gpV)
    print gV*gV, "-", gW*gW

    X, Y, F = gV.export(0.0, 90.0, 50, 0.0, 40.0, 50)
    X, Y, G = gW.export(0.0, 90.0, 50, 0.0, 40.0, 50)
    if i==0:
        f = plt.figure(figsize=(14, 6), dpi=120)
    else:
        f = plt.gcf()
    f.clear()
    ax = f.add_subplot(121, projection='3d')
    ax.plot_wireframe(x, y, abs(F))
    bx = f.add_subplot(122, projection='3d')
    bx.plot_wireframe(x, y, abs(G))
    plt.tight_layout()
    plt.draw()
    #plt.show()
t1 = time.time()

print t1 - t0

plt.ioff()
plt.show()
