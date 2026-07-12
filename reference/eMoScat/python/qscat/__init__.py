# First append the location of C++ library to the python path

import sys, os
sys.path.append("/opt/intel/composer/lib/intel64")
import ctypes
try:
	ctypes.CDLL('libmkl_rt.so', ctypes.RTLD_GLOBAL)
except OSError:
	print('Intel module not found, continuing...')


fdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(fdir)
from pyArrays import *
