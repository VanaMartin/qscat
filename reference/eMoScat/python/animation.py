#!/usr/bin/env python3

"""
    Python Animation saver
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#from pylab import *
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

out = np.zeros(E.shape, np.complex) + 1e-90

ierg = -0.0976049
#esh = ierg  + 0.125 #-0.744777
dt = T[1] - T[0]
delta = 0 if not options.elastic else -1.0

S = np.zeros(E.shape, np.complex) + 1e-90

#for i in range(0, E.shape[0]):
#    if abs(F[i]) > 0:
#        e = E[i] + ierg;
#        eit = np.exp(1j * e * T )
#        val = np.inner(eit, C)
#        val  = val * dt / (F[i] * F0[i]) / ( 2.0 * pi)
#        out[i] = abs((val-delta)/ (2*pi))**2 * 4 * pi**3 / (2 * E[i])
#    else:
#        out[i] = np.nan

dpi = 200

#def update_img(n):
#    tmp = rand(300,300)
#    im.set_data(tmp)
#    return im

W = dt / (F * F0) / ( 2.0 * pi)

print(plt.get_backend())
plt.switch_backend('TKAgg')
matplotlib.rc('text',usetex=True)

def ani_frame(T,C,W,S,erg):

    plt.ion()
    fig = plt.figure()
    fig.suptitle("H$_2^+$+ $\Sigma$ symmetry dissociative recombination cross sections VE 0 $\\rightarrow$ DR 1")
    ax = fig.add_subplot(111)

    ln, = ax.plot(E, np.real(out))
    fig.set_size_inches([10,5])
    ax.set_ylim([1e-14,1e+2])
    ax.set_xlim([0,0.05])
    ax.set_yscale('log')
    #plt.show()
    #fig.tight_layout()
    #fig.canvas.draw()
    ax.set_title('T = 0')


    def update_img(n):
        global S
        if n%100 == 0:
            print(n)
        k = 10
        for i in range(0,k):
            ceit = C[n*k+i] * np.exp( (E + erg) * 1j * T[n*k+i] )
            S += W * ceit
        out = abs((S-delta)/ (2*pi))**2 * 4 * pi**3 / (2 * E)
        ln.set_data(E, np.real(out))
        ax.set_title('T = %d' % (n * k * int(dt)) )
        return ln


    ani = animation.FuncAnimation(fig,update_img,50000, interval=20)
    #writer = animation.writers['ffmpeg'](fps=30)

    #for i in range(0,100000):
    #    update_img(i)
    #    fig.canvas.draw()

    ani.save('h2p_10x.mp4', writer = 'avconv',dpi=dpi)
    #return ani
    plt.ioff()
    return

ani_frame(T,C,W,S,ierg)
