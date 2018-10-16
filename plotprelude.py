import matplotlib
import numpy
from matplotlib import pylab

pylab.rc("savefig", dpi=300)
pylab.rc("path", simplify=True)
pylab.rc("font", family="serif")
pylab.rc("mathtext", fontset ="custom")
pylab.rc("text", usetex=False)


for x in ["xtick", "ytick"]:
    pylab.rc(x+".major", size=6, width=1)
    pylab.rc(x+".minor", size=3, width=1)

pylab.rc("lines", markeredgewidth= 1)

pylab.rc("legend", numpoints=1, frameon= False, handletextpad=   0.3)

plotstyle={ "capsize": 4, "capthick": 0.5, "elinewidth": 0.5}


matplotlib.style.use('seaborn-paper')

def mkfig():
    return pylab.figure(figsize=(5., 4.), dpi=200)

def pp_ps(d, title=""):
    fig=mkfig()
    ax=fig.add_subplot(111)
    xvar = numpy.arange(d.shape[0])-(d.shape[0]+1)/2
    pylab.semilogy(xvar,d)
    ax.set_xlabel("Delay ( %f ns)" % (1e9/(100e6/1024) / 128,))
    ax.set_ylabel("Absolute mean cross-spectral power density ")
    ax.set_title("PowerSpectrum "+title )
    ax.xaxis.set_ticks_position('bottom')
    ax.grid()
