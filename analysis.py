import numpy


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
