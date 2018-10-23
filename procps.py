"""
Power Spectrum Estiation from closure phase data

General notes:

dc: Data Cube of phases of dimensions:
    time x day x triad x freq 
"""

import numpy
import scipy, scipy.signal


def cmedn(a, axis, n):
    """Circular MEDian N

    Returns n values around the median computed along axis using
    circular stats. NB values are not necessarily sorted in distance
    from median (argpartition is used)
    """
    p=a[numpy.newaxis, :]-a[:, numpy.newaxis]
    p=(p+numpy.pi) % (2 * numpy.pi) - numpy.pi
    p = numpy.abs(p).sum(1)
    ii=numpy.argpartition(p, n+1, axis=axis)
    return numpy.choose(ii[0:n], a)

def mdays(dc, n):
    """Median across Days

    Compute cmedn across the days axis. Needed because cmedn is N^2
    memory scaling
    """
    dnew=numpy.array([[cmedn(dc[i,:,j], axis=0, n=n) for j in range(dc.shape[2])] for i in range(dc.shape[0])])
    # Rotate axes so days are in position 1 again
    return numpy.moveaxis(dnew, 2, 1)

def psX(a, b):
    return scipy.signal.csd(a, b,
                            nperseg=128,
                            detrend=None)

def psXtime(dcm, ii):
    return psX(dcm[ii:,0],
               dcm[:-ii,0])

def psXtimeCTimeITri(dcm, ii):
    ps=psXtime(dcm, ii)
    CTime=numpy.mean(ps[1],
                     axis=0)


    ITri=numpy.mean(numpy.abs(CTime),
                    axis=0)
    return ITri

def psXtimeCTimeCTriR(dcm, ii):
    ps=psXtime(dcm, ii)
    CTime=numpy.mean(ps[1],
                     axis=0)
    CTri=numpy.mean(CTime,
                    axis=0)
    return CTri.real


def psXmedCTimITri(dcm):
    ps=scipy.signal.csd(dcm[:,0],
                        dcm[:,1],
                        nperseg=128,
                        detrend=None)    
    CTim=numpy.mean(ps[1], axis=0)
    ITri=numpy.mean(numpy.abs(CTim),
                    axis=0)
    return ITri

def psXmedCTimCTri(dcm, window="hann"):
    ps=scipy.signal.csd(dcm[:,0],
                        dcm[:,1],
                        nperseg=128,
                        detrend=None,
                        window=window)    
    CTim=numpy.mean(ps[1], axis=0)
    CTri=numpy.mean(CTim,
                    axis=0)
    return numpy.abs(CTri)

def psXmedCTimCTriR(dcm, window="hann"):
    ps=scipy.signal.csd(dcm[:,0],
                        dcm[:,1],
                        nperseg=128,
                        detrend=None,
                        window=window)    
    CTim=numpy.mean(ps[1], axis=0)
    CTri=numpy.mean(CTim,
                    axis=0)
    return CTri.real

def psXTrCTimCTr(dcm, window="hann",
                 shuffle=True):
    dcm=dcm.copy()
    if shuffle:
        #Scramble along the triad axis
        dcm=shuffleax(dcm, 2)
    ps=scipy.signal.csd(dcm[:,0,1:],
                        dcm[:,0,:-1],
                        nperseg=128,
                        detrend=None,
                        window=window)    
    CTim=numpy.mean(ps[1], axis=0)
    CTri=numpy.mean(CTim,
                    axis=0)
    return CTri

def lmodel(fin):
    """Load model as produced by cclosure"""
    d=numpy.load(fin)
    cc=d["phase"]
    # Rotate axes so that the shape is (time x 1 x triad x channel)
    cc=numpy.moveaxis(numpy.moveaxis(cc, 1, 0), -1, 0)
    return cc

def shuffleax(d, ax):
    """Shuffle along axis ax"""
    dd=numpy.moveaxis(d, ax, 0)
    numpy.random.shuffle(dd)
    return numpy.moveaxis(dd, 0, ax)
