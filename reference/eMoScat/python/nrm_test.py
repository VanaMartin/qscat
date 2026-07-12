#!/usr/bin/env python3

 # First append the location of C++ library to the python path
import sys
sys.path.append("/opt/intel/composer/lib/intel64")
import ctypes
try:
	ctypes.CDLL('libmkl_rt.so', ctypes.RTLD_GLOBAL)
except OSError:
	print('Intel module not found, continuing...')

import qscat as qs
import numpy as np

import json, math, cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d

def expandXY(X,Y):
    XX = np.zeros([X.shape[0],Y.shape[0]])
    YY = np.zeros([X.shape[0],Y.shape[0]])

    for i in range(0,Y.shape[0]):
        XX[:,i] = X

    for i in range(0,X.shape[0]):
        YY[i,:] = Y
    return XX, YY

def check(E,i,idx, c=0):
    q = abs(E[idx,i] - E[idx,i+1]) #/E[idx, i+1]
    r = abs(E[idx-1,i] - E[idx,i+1])
    if q > r:
        print(q, c)
        idx += 1
        c += 1
        return check(E,i,idx,c)
    return idx

def CPT(E,i,idx):
    idx = check(E,i,idx)

    if idx > 1:
        k = idx - 1
        R = np.zeros([k],np.complex)
        R[:] = E[0:k, i]
        E[0,i] = E[k,i]
        E[1:idx,i] = R[:]

    return E, idx


with open("input/experimental/N2-model.json") as f:
    try:
        c = "".join(f.readlines())
        p = json.loads(c)
    finally:
        f.close()

gx = qs.gridFromString(json.dumps(p["grids"]["electronic_short"]))
gy = qs.gridFromString(json.dumps(p["grids"]["nuclear"]))

g = qs.grid2d(gx,gy)

phi = qs.gVector2d(g)
aux = qs.gVector(gx)
#pot = qs.gVector2dPotential(g, json.dumps(p["model"]["potential"]))

pstr = json.dumps(p["model"]["potential"])

V0 = qs.gVectorMorse(gy, pstr)

E = np.zeros([gx.nb(), gy.nb()], np.complex)

y = np.real(gy.export())
x = np.real(gx.export())

smooth = 1.0 - 1.0 / (1.0 + np.exp(-x+10))

print(gy.nb())
idx = 1

print("Loading phi_d")
if not phi.load("phid.qbin"):
    print("Not successful!")
    for i in reversed(range(0, gy.nb())):
        p = qs.gVectorElectronicLambda(gx, gy.x(i), pstr)
        H = qs.fOperator(gx, gx.nb(), 0).addKineticTerm(1.0)
        H += p
        eSys = H.eigen_system()
        #print(eSys.eigenvalue(0))
        E[:,i], psi = eSys.export()

        if i < gy.nb()-2:
            E, idx = CPT(E,i,idx)

        b = eSys.state(idx-1)
        for j in range(0,gx.nb()):
            b[j] *= smooth[j]
        aux = qs.gVector(gx, b, False)
        nrm = math.sqrt(abs(aux*aux))
        alpha = cmath.phase(aux[20])
        aux *= 1.0/nrm * cmath.exp( -1j * alpha )
        phi.write_x(aux, i)


phi.save("phid.qbin")

if False:
    plt.plot(y[0:200], E[:,0:200].transpose(), marker='x')
    plt.show()


if True:
    f = plt.figure()
    ax = f.add_subplot(111, projection='3d')

    Y,X,Z = phi.export(0, 20, 200, 0.5, 4, 200)

    XX, YY = expandXY(X,Y)

    ax.plot_surface(XX,YY,np.abs(Z)**2, rcount=100, ccount=100)
    plt.show()

exit(0)

Vd = qs.gVector(gy)

Vdk = [ qs.gVector(gx) for i in range(0,gx.nb()) ]
idx = range(0,gx.nb())

for i in reversed(range(0, gy.nb())):
    print("R =", gy.xr(i), "...")
    Hel = qs.fOperator(gx)
    Hel.addKineticTerm(1.0)
    Vel = qs.gVectorElectronicLambda(gx, gy.x(i), pstr)
    Hel += Vel

    p = phi.read_x(i)
    Q = qs.fOperator(gx)
    Q.outer_product(p,p)

    P = qs.fOperator(gx)
    P.id()
    P -= Q

    PHP = P.copy()
    PHP *= Hel
    PHP *= P

    a = qs.gVector(gx)
    PHP.gemv(1.0, p, 0.0, a)

    val = V0.getVal(i) + p*a
    print(val, p*a)
    Vd.setVal( val, i )

    aH = PHP.copy()
    eSys = aH.eigen_system()

    E[:,i], psi = eSys.export()
    if i < gy.nb() - 3 and True:
        for j in reversed(range(1,gx.nb())):
            h = (gy.x(i+1) - gy.x(i+2)).real
            a = (E[j,i+1] - E[j,i+2]).real
            g = h *(gy.x(i)-gy.x(i+1)).real
            e = E[j,i+1].real + a / h * g
            #if abs(E[j,i] - e) > abs(1.0e-2 * g):
            R = np.zeros([j],np.complex)
            R[:] = E[0:j,i].real - e
            jj = np.argmin(R)
            if j != jj:
                print("...", abs(E[j,i] - e), R[jj], ":", jj, "<->", j)
                e0 = E[j,i]
                E[j,i] = E[jj,i]
                E[jj,i] = e0

if True:
    f = plt.figure()
    ax1 = f.add_subplot(121)
    ax1.plot(y[0:300], np.real(E[:,0:300]).transpose(), marker='x')
    ax2 = f.add_subplot(122, sharex=ax1)
    ax2.plot(y[0:300], np.imag(E[:,0:300]).transpose(), marker='x')
    plt.show()

X,Y = Vd.export(0,10, 200)

if False:
    plt.plot(X,Y)
    plt.show()
