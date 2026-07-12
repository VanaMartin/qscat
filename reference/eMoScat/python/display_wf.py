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
import matplotlib.pyplot as plt
from optparse import OptionParser

""" ENTERING EXPERIMENTAL ZONE

This is an attempt to utilize matplotlib colors and image modules to build a
fully functional 2D color mapping. The main purpose of this module is to
provide an easy to use method for displaying complex values in two dimension.
However the extensions to some other cases (such as normal vectors) are
currently in consideration.  """

import math
import matplotlib
from matplotlib import cm
import matplotlib.colors as mcolors

# Basic utility for converting the complex numpy array to HSV model
def complex_to_hsv(c, Mag=1.0, inverse=False):
    """ The complex ndarray to HSV model conversion. Note that the sue
    conversion seems unnecessairly complicated, however there is a good reason.
    The hue value should be zero for zero angle to maintain the consistency
    with other implementations."""
    # first ensure np.array
    c = np.asarray(c)

    hsv = np.zeros(c.shape + (3,))

    h = hsv[..., 0]
    s = hsv[..., 1]
    v = hsv[..., 2]

    if not np.iscomplexobj(c):
        raise ValueError("The input array must be of compex number "
                         "dtype {dtp} was found.".format(dtp=c.dtype))

    h[:] = (np.angle(c)) / (2*math.pi)                                      # complex angle, but on [-pi,pi)
    h[:] = np.piecewise(h, [h<0,h>0], [lambda x: 1.0 + x, lambda x: x] )    # solve the negative part via piecewise function
    if inverse:
        s[:] = np.abs(c) / Mag                                                  # absolute value divided by magintude
        v[:] = 1.0 / np.maximum(1.0, s)                                         # saturation (1.0 if abs(c) < Mag)
        s[:] = np.minimum(1.0, s)                                               # value (1.0 if abs(c) > Mag)
    else:
        v[:] = np.abs(c) / Mag                                                  # absolute value divided by magintude
        s[:] = 1.0 / np.maximum(1.0, v)                                         # saturation (1.0 if abs(c) < Mag)
        v[:] = np.minimum(1.0, v)                                               # value (1.0 if abs(c) > Mag)

    # we are done, returning hsv array
    return hsv

def complex_to_rgb(c, Mag=1.0):
    return mcolors.hsv_to_rgb(complex_to_hsv(c,Mag))

def complex_to_irgb(c, Mag=1.0):
    return mcolors.hsv_to_rgb(complex_to_hsv(c,Mag,inverse=True))

def complex_to_rgba(c, Mag=1.0):                # FIXME
    return mcolors.to_rgba(complex_to_rgb(c,Mag))

############ MAIN ###############3

title = "Indirect scattering of H${_2}^{+}$ with electron "
pot_name = "H2+ potential in sigma symmetry"

matplotlib.rc('text',usetex=True)
### PARSING ARGUMENTS
parser = OptionParser(usage="%prog [options] input/dir")
parser.add_option("-i", "--inverse", action='store_true', dest="inv", help="color inversion")
parser.add_option("-m", "--magnitude", action='store', dest="mag", help="sets the HSL scaling magnitude", metavar="float")
parser.add_option("-f", "--zoom-factor", action='store', dest="zoomf", help="Zooming factor for region zoom", metavar="float")

(options, args) = parser.parse_args()

x = [0, 50]
y = [0, 10]

sampling = 800
scaling = 4.5
factor = 10.0 if not options.zoomf else float(options.zoomf)
mag = 0.005 if not options.mag else float(options.mag)

fname = args[0]
inverse = options.inv if options.inv else False

psi = pa.gVector2d()
if not psi.load(fname):
    print("Cannot read file")
    exit(-1)

pot = pa.gVector2d()
if not pot.load("pot.qbin"):
    print("Cannot read file")
    exit(-1)

X,Y,P = pot.export(x[0],x[1],sampling, y[0],y[1],sampling)
X,Y,Z = psi.export(x[0],x[1],sampling, y[0],y[1],sampling)

border=X.shape[0] - 1
for i in range(0,X.shape[0]):
    if X[i] >= scaling:
        border = i-1
        break

Z[:,border:] *= factor

#plevels = np.arange(-10,10,0.1)
pp = np.arange(0,2,0.05)
plevels = np.exp(pp) - 1
plevels = np.append(-plevels[::-1],plevels[1:])

wlevels = np.arange(0.0, 4*mag, mag/5)


if inverse:
    ctcolor = "k"
    rgba = complex_to_irgb(Z, mag)
else:
    ctcolor = "w"
    rgba = complex_to_rgb(Z, mag)

plt.ion()
f = plt.figure()
ax = f.add_subplot(111)
f.tight_layout()
ax.axes.imshow(rgba, origin='upper', interpolation='bessel', aspect='auto', extent=[y[0],y[1],x[1],x[0]])
matplotlib.rcParams['contour.negative_linestyle']='solid'
pcmap = None
pc = ax.contour(X,Y,np.real(P), levels=plevels, colors='grey', alpha=0.5)
plt.setp(pc.collections[int(plevels.shape[0]/2)], linestyle="dashed", color=ctcolor)

ax.contour(X,Y,np.abs(Z)**2, levels=wlevels**2, colors=ctcolor, alpha=0.5)
ax.plot([X[border],X[border]],[Y[0],Y[-1]], linewidth=5.0, color=ctcolor)

plt.xlabel("Nuclear coordinate R [a$_0$]")
plt.ylabel("Electronic coordinate r [a$_0$]")
plt.title(title)

#ax.xaxis.set_color(ctcolor)
ax.tick_params(color=ctcolor)

### Legend generator
## LEGEND
sx = np.arange(-2,2,0.1) * mag
smp = np.add.outer(sx*1.0j, sx)
if inverse:
    lrgb = complex_to_irgb(smp,mag)
else:
    lrgb = complex_to_rgb(smp,mag)

sax = plt.axes([.7, .2, .2, .2], axisbg='y')
sax.axes.imshow(lrgb, origin='upper', interpolation='bessel', extent=[-2*mag,2*mag,-2*mag,2*mag])
sax.contour(sx,sx,np.abs(smp)**2, levels=wlevels**2, colors=ctcolor, alpha=0.5)
sax.set_xlabel("Re(z)")
sax.set_ylabel("Im(z)")

tcolor = ctcolor
sax.spines['top'].set_color(tcolor)
sax.spines['bottom'].set_color(tcolor)
sax.spines['left'].set_color(tcolor)
sax.spines['right'].set_color(tcolor)
sax.tick_params(color=tcolor)
sax.xaxis.label.set_color(tcolor)
sax.yaxis.label.set_color(tcolor)
for t in sax.xaxis.get_ticklabels():
    t.set_color(tcolor)
for t in sax.yaxis.get_ticklabels():
    t.set_color(tcolor)

ax.text(5.5, 35, "H$_2^+$ $\Sigma$ symmetry", fontsize=22, weight='bold', color=tcolor, zorder=10)

ax.text(8.5, 4.75, "Dissociative recombination", fontsize=16, weight='bold', color=tcolor, zorder=10)
ax.text(2.3, 48, "Vibrational excitation", fontsize=16, weight='bold', color=tcolor, zorder=10)
ax.arrow(8.5, 3.0, 1.0, 0, linewidth=3.0, head_width=0.8, head_length=0.2, fc=tcolor, ec=tcolor, zorder=10)
ax.arrow(2.0, 42, 0, 5, linewidth=3.0, head_width=0.08, head_length=2.0, fc=tcolor, ec=tcolor, zorder=10)

plt.ioff()
plt.show()
