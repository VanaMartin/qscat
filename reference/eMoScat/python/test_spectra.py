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

with open("input/experimental/N2-model.json") as f:
    try:
        c = "".join(f.readlines())
        p = json.loads(c)
    finally:
        f.close()

gx = qs.gridFromString(json.dumps(p["grids"]["electronic"]))
gr = qs.gridFromString(json.dumps(p["grids"]["electronic_short"]))

r = gr.export()
x = gx.export()

plt.plot(np.real(x), np.imag(x))
plt.plot(np.real(r), np.imag(r))
plt.show()

y = 1.5

pstr = json.dumps(p["model"]["potential"])
px = qs.gVectorElectronicLambda(gx, y, pstr)
pr = qs.gVectorElectronicLambda(gr, y, pstr)

Hx = qs.fOperator(gx, gx.nb(), 0).addKineticTerm(1.0)
Hx += px
Hr = qs.fOperator(gr, gr.nb(), 0).addKineticTerm(1.0)
Hr += pr

Ex = Hx.eigen_system()
Er = Hr.eigen_system()

ex, psix = Ex.export()
er, psir = Er.export()

f = plt.figure()
ax = f.add_subplot(111)

ax.plot(np.real(ex), np.imag(ex), marker="o", linewidth=0.0)
ax.plot(np.real(er), np.imag(er), marker="x", linewidth=0.0)

ax.set_xlim([0,1.5])
ax.set_ylim([-0.8,0.05])

plt.show()

print("got it")
exit(0)
