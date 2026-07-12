#!/usr/bin/env python3

 # First append the location of C++ library to the python path
import sys
sys.path.append("/opt/intel/composer/lib/intel64")
import ctypes
try:
	ctypes.CDLL('libmkl_rt.so', ctypes.RTLD_GLOBAL)
except OSError:
	print('Intel module not found, continuing...')

# Python Arrays Wrapper

import pyArrays as arr
import numpy as np
import math

# Vector <int>
print("Integer")
A = arr.iVector(4)
B = arr.iVector(4)
A.fill(1)
B.fill(5)

print(B[0])
B[0] = 10

print(B[0])

C = arr.iVector(B)
B += A

print(B.export(), C.export())

# vector <double>
print("Double case")
A = arr.dVector(4)
A.fill(2.0)

print(A.export())

B =  A * 3.0

print(B.size())
print(A*B, A.export(), B.export())

# vector <complex>
print("Complex case")

Z = arr.zVector(10)

Z.fill(1.0+1.0j)

print(Z.export(), Z*Z)

# matrix <double>
print("Double matrix")
Q = arr.dMatrix(4,4).fill(0.5)
print(Q.export())

Id = arr.dMatrix().set_identity(4)
print(Id.export())
Q += Id

A.fill(2.0)
L =  Q*A

print(A.export(), L.export())

L = Q.solve(A)

print(A.export(), L.export())

print("Complex matrix")
Z = arr.zMatrix().set_identity(4)
for i in range(0,4):
	if i < 3 : Z[i*5 + 1] = -0.5
	if i > 0 : Z[i*5 - 1] = -0.5

print(Z.export())

print("Eigensystem")
E, W = Z.eigen_system().export()

print(E)
print(W)
