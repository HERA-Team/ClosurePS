import numpy

if 0:
    finXX=numpy.load("/bigfast/temp1/hera2/jck/EQ14binnedXX.npz")
    XXm=mdays(finXX["closures"], 1)

if 0:
    pp_wfl(finXX["closures"][:,0,0], finXX["last"][:,0],
           title="Day 0, triad 0")
    pylab.savefig("plots/wfl_day0.png", bbox_inches='tight')
    
    pp_wfl(XXm[:,0,0], finXX["last"][:,0],
           title="Median filtered across all days, triad 0")
    pylab.savefig("plots/wfl_whole_med.png", bbox_inches='tight')

    pp_wfl(XXm[:,0,0,500:750], finXX["last"][:,0],
           title="Median filtered across all days, triad 0",
           vrange=(-0.5,0.5),
           gg=True,
           choff=500)
    pylab.savefig("plots/wfl_ribs.png", bbox_inches='tight')

    XXms=XXm[150:250, :, 0, 120:380]

if 0:
    pp11=mdays(finXX["closures"][:,0:9], 1)
    pp12=mdays(finXX["closures"][:,9:], 1)
    pp33=numpy.concatenate((pp11, pp12), axis=1)

if 0:
    pylab.clf()
    pp_cps([ (pp11[208, 0, 0, 120:900], "First half of days"),
             (pp12[208, 0, 0, 120:900], "Second half of days")],
           choff=120)
    pylab.ylim((-0.75, 0.75))
    pylab.savefig("plots/specs_median_halves.png")

if 0:
    finEWXX=numpy.load("/bigfast/temp1/hera2/jck/LinEWXX.npz")
    EWXXm=mdays(finEWXX["closures"],1)

if 1:
    for ii in [2,4,6,8]:
        pylab.clf()
        x=EWXXm[175:(175+4*ii),0,0,120:900]
        pp_cps([ (x.reshape( (-1, ii, 780)).mean(axis=1)[i], "record %i" % i) for i in range(4)],
               choff=120,
               title="%i consecutive records average" % ii)
        pylab.ylim((-0.75, 0.75))
        pylab.savefig("plots/specs_tstep-%i.png" % ii , bbox_inches='tight')

if 0:
    fin=numpy.load("/bigfast/temp1/hera2/jck/EQ14binnedXX.npz")
    dd=fin["closures"]
    pp2=mdays(dd,2)

if 0:
    pp11=mdays(dd[:,0:8], 1)
    pp12=mdays(dd[:,8:], 1)
    pp33=numpy.concatenate((pp11, pp12), axis=1)

if 0:
    pp41=mdays(dd[:,::2], 1)
    pp42=mdays(dd[:,1::2], 1)
    pp4=numpy.concatenate((pp41, pp42), axis=1)    
    
if 0:
    x1=psXmedCTimCTriR(numpy.exp(1j*(pp2[150:151,:,10:11,120:380])),
                       window="hamming")
    pp_ps(x1)
    pylab.savefig("plots/spec-noisedom.png")

if 0:
    x1=psXmedCTimCTriR(numpy.exp(1j*(pp33[150:250,:,5:6,120:380])),
                       window="hamming")
    pp_ps(x1)
    pylab.savefig("plots/spec-onetriad-twomedian.png")    

if 0:
    x1=psXmedCTimCTriR(numpy.exp(1j*(dd[150:250,:,:,120:380])),
                       window="hamming")
    pp_ps(x1)
    pylab.savefig("plots/spec-noisedom-nomed.png")    

if 0:
#Track down the coherence
    pylab.clf()
    pylab.plot(pp2[150,0,10,120:380])
    pylab.plot(pp2[150,1,10,120:380])
    pylab.savefig("plots/whatgoesin.png")

    pylab.clf()
    pylab.plot(dd[150,0,10,120:380])
    pylab.plot(dd[150,1,10,120:380])
    pylab.savefig("plots/whatgoesin-nomed.png")

    pylab.clf()
    pylab.plot(pp33[150,0,10,120:380])
    pylab.plot(pp33[150,1,10,120:380])
    pylab.savefig("plots/whatgoesin-twomedian.png")

if 0:    
    pylab.clf()
    pylab.plot(pp4[150,0,10,120:380])
    pylab.plot(pp4[150,1,10,120:380])
    pylab.savefig("plots/whatgoesin-inteleavedmedian.png")            

if 0:
    pylab.clf()
    pylab.plot(pp33[150,0,10,120:380])
    pylab.plot(pp33[151,0,10,120:380])
    pylab.savefig("plots/whatgoesin-tworecs.png")
    

if 0:

    FR=550,750
    pylab.clf()
    pylab.plot(pp33[200,0,10,FR[0]:FR[1]])
    pylab.plot(pp33[200,1,10,FR[0]:FR[1]])
    pylab.savefig("plots/whatgoesin-twomedian-hf.png")

    pylab.clf()
    pylab.plot(pp33[200,0,10,FR[0]:FR[1]])
    pylab.plot(pp33[201,0,10,FR[0]:FR[1]])
    pylab.savefig("plots/whatgoesin-tworecs-hf.png")

# Look for x-corr between triads
if 0:
    fin=numpy.load("/bigfast/temp1/hera2/jck/EQ14binnedYY.npz")
    dd=fin["closures"]
    pp11=mdays(dd[:,0:8], 1)
    pp12=mdays(dd[:,8:], 1)
    pp33=numpy.concatenate((pp11, pp12), axis=1)

if 0:
    dw=numpy.exp(1j*pp33[150:250, :, :, 120:380])

if 0:    
    dwp=psX(dw[:, :, 0, 120:380],
            dw[:, :, 20, 120:380])
    pp_ps(numpy.transpose(dwp[1][:,0]))    
    pylab.savefig("plots/trxcor-1.png")

    p1=numpy.mean(dwp[1], axis=0)
    print (p1.shape)
    pp_ps(p1[0])    
    pylab.savefig("plots/trxcor-2.png")        


if 0:
    pylab.clf()
    pylab.plot(dw[0,0,0])
    pylab.plot(dw[0,0,1])
    pylab.savefig("plots/tr01.png")

    pylab.clf()
    pylab.plot(dw[0,0,0])
    pylab.plot(dw[0,0,10])
    pylab.savefig("plots/tr10.png")


if 0:
    pp_ps([(numpy.abs(psXTrCTimCTr(dw)), "scrambled"),
           (numpy.abs(psXTrCTimCTr(dw, shuffle=False)), "ordered")])
    pylab.savefig("plots/psXTrCTimCTr.png")

    x1=psXmedCTimCTri(dw)
    pp_ps(x1)
    pylab.savefig("plots/psXmedCTimCTri.png")    


if 0:
    fin=numpy.load("/bigfast/temp1/hera2/jck/EQ14binnedYY.npz")
    dd=fin["closures"]
    pp11=mdays(dd[:,0:8], 1)
    pp12=mdays(dd[:,8:], 1)
    shuffleax(pp12, 2)    
    pp33=numpy.concatenate((pp11, pp12), axis=1)

if 0:
    dw2=numpy.exp(1j*pp33[150:250, :, :, 120:380])

    # cross triand cross days
    x2=psXmedCTimCTri(dw2,
                      window="hamming")
    pp_ps([ (numpy.abs(x1), "first half X second half same triad"),
            (numpy.abs(x2), "first half X second half random  triad")])
    pylab.savefig("plots/psec-xdayxtriad.png")        

    
    
