import matplotlib
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
    return pylab.figure(figsize=(5., 4.), dpi=300)

def pp_ps(d, pf=""):
    fig=mkfig()
    ax=fig.add_subplot(111)
    if type(d) is list:
        for dd,dl in d:
            pylab.semilogy(dd, label=dl)
            ax.legend()
    else:
        pylab.semilogy(d)
    ax.set_xlabel("Delay ( %f ns)" % (1e9/(100e6/1024) / 128,))
    ax.set_ylabel("Absolute mean cross-spectral power density ")
    ax.set_title("PowerSpectrum  " )
    ax.xaxis.set_ticks_position('bottom')
    ax.grid()

def pp_wfl(dc, lst,
           title="",
           vrange=(None, None),
           gg=False,
           choff=0):
    fig=pylab.figure(figsize=(8., 4.), dpi=300)
    ax=fig.add_subplot(111)    
    mappbl=ax.matshow(dc,
                      vmin=vrange[0],
                      vmax=vrange[1])
    ax.set_xlabel("Channel Number")
    ax.set_ylabel("LST (hours from midnight) ")
    ax.set_title("Closure phase waterfall -- %s" % title)
    yticks = ["{:2.2f}".format((l %1 ) * 24) for l in numpy.modf(lst)[0]]
    ax.set_yticks(numpy.arange(len(yticks))[::50] )
    ax.set_yticklabels(yticks[::50])

    xticks = numpy.arange(dc.shape[1])[::50]
    ax.set_xticks(xticks)
    ax.set_xticklabels(["%i" % (xx+choff) for xx in xticks])  
    ax.xaxis.set_ticks_position('bottom')
    if gg:
        ax.grid()    
    cb=pylab.colorbar(mappbl)
    cb.set_label("Closure phase (rad)")
    

    
