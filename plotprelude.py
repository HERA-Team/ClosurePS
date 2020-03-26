import matplotlib
from matplotlib import pyplot

pyplot.rc("savefig", dpi=300)
pyplot.rc("path", simplify=True)
pyplot.rc("font", family="serif")
pyplot.rc("mathtext", fontset ="custom")
pyplot.rc("text", usetex=False)


for x in ["xtick", "ytick"]:
    pyplot.rc(x+".major", size=6, width=1)
    pyplot.rc(x+".minor", size=3, width=1)

pyplot.rc("lines", markeredgewidth= 1)

pyplot.rc("legend", numpoints=1, frameon= False, handletextpad=   0.3)

plotstyle={ "capsize": 4, "capthick": 0.5, "elinewidth": 0.5}


matplotlib.style.use('seaborn-paper')

def mkfig():
    return pyplot.figure(figsize=(5., 4.), dpi=300)

def pp_ps(d, pf=""):
    fig=mkfig()
    ax=fig.add_subplot(111)
    if type(d) is list:
        for dd,dl in d:
            pyplot.semilogy(dd, label=dl)
            ax.legend()
    else:
        pyplot.semilogy(d)
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
    fig=pyplot.figure(figsize=(8., 4.), dpi=300)
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
    cb=pyplot.colorbar(mappbl)
    cb.set_label("Closure phase (rad)")
    
def pp_cps(dl,
           choff=0,
           title="",
           alpha=0.5):
    fig=pyplot.figure(figsize=(8., 4.), dpi=300)
    ax=fig.add_subplot(111)
    if type(dl) is list:
        for dd,dll in dl:
            pyplot.plot(dd, label=dll,
                       marker=".",
                       linestyle="-",
                       linewidth=0.5,
                       alpha=alpha)
        ax.legend()
        xticks = numpy.arange(dl[0][0].shape[0])[::100]
    ax.set_xticks(xticks)
    ax.set_xticklabels(["%i" % (xx+choff) for xx in xticks])  
    ax.xaxis.set_ticks_position('bottom')        
    ax.set_xlabel("Channel Number")
    ax.set_ylabel("Closure phase (rad)")
    ax.set_title("Closure phase spectrum -- %s" % title)    

def pp_avar(ddl, fnameout, tl,
            inttime=10.7668):
    lines = ["-","--","-.",":"]
    linecycler = cycle(lines)
    pyplot.close("all")
    f, axl=pyplot.subplots(len(ddl),
                          sharex=True,
                          sharey=False,
                          figsize=(8, 4), dpi=300)
    if len(ddl)==1: axl=[axl]
    for i,dd in enumerate(ddl):
        ax=axl[i]
        for d,l in dd:
            line,=ax.plot(d[0]*inttime, numpy.rad2deg(d[1]), next(linecycler), label=l)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_title(tl[i])
    ax.set_xlabel("Stride and integration time (s)")
    axl[len(axl)//2].set_ylabel("Allan deviation (degrees closure phase)")
    #ax.set_xticks(d[0]*inttime)
    for tax in [ax.get_xaxis(), ax.get_yaxis()]:
        tax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    f.legend(fancybox=True)
    f.savefig(fnameout, bbox_inches='tight')
    
