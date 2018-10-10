import numpy
import os
import glob
import casac

n_channels = 1024

mslist = glob.glob("*ms")
mslist = sorted(mslist)
print(mslist)


#### List of Baselines I care about

ew_14m = [[0,1],[1,2],[11,12],[12,13],[13,14],[23,24],[24,25],[25,26],[26,27],[36,37],[37,38],[38,39],[39,40],[40,41],[51,52],[52,53],[53,54],[54,55],[65,66],[66,67],[67,68],[68,69],[69,70],[70,71],[82,83],[83,84],[84,85],[85,86],[86,87],[87,88],[120,121],[121,122],[122,123],[123,124]]

ew_28m = [[0,2],[11,13],[12,14],[23,25],[24,26],[25,27],[36,38],[37,39],[38,40],[39,41],[51,53],[52,54],[53,55],[65,67],[66,68],[67,69],[68,70],[69,71],[82,84],[83,85],[84,86],[85,87],[86,88],[120,122],[121,123],[122,124]]

# Origin on RHS of array
deg60_14m = [[0,12],[1,13],[2,14],[11,24],[12,25],[13,26],[14,27],[23,37],[24,38],[25,39],[26,40],[27,41],[36,51],[37,52],[38,53],[39,54],[40,55],[51,67],[52,68],[53,69],[54,70],[55,71],[65,82],[66,83],[67,84],[68,85],[69,86],[70,87],[71,88],[120,140],[121,141],[122,142],[123,143]]

deg60_28m = [[0,25],[1,26],[2,27],[11,38],[12,39],[13,40],[14,41],[23,52],[24,53],[25,54],[26,55],[36,67],[37,68],[38,69],[39,70],[40,71],[51,84],[52,85],[53,86],[54,87],[55,88]]

deg120_14m = [[0,11],[1,12],[2,13],[11,23],[12,24],[13,25],[14,26],[23,36],[24,37],[25,38],[26,39],[27,40],[37,51],[38,52],[39,53],[40,54],[41,55],[51,66],[52,67],[53,68],[54,69],[55,70],[66,82],[67,83],[68,84],[69,85],[70,86],[71,87],[120,139],[121,140],[122,141],[123,142],[124,143]]

deg120_28m = [[0,23],[1,24],[2,25],[11,36],[12,37],[13,38],[14,39],[24,51],[25,52],[26,53],[27,54],[36,65],[37,66],[38,67],[39,68],[40,69],[41,70],[51,82],[52,83],[53,84],[54,85],[55,86]]

listoflistoftriads = []


listoflistoftriads.append(ew_14m)
listoflistoftriads.append(ew_28m)
listoflistoftriads.append(deg60_14m)
listoflistoftriads.append(deg60_28m)
listoflistoftriads.append(deg120_14m)
listoflistoftriads.append(deg120_28m)



snrmap = numpy.zeros(shape=(6,1024,60*len(mslist)),dtype=numpy.complex64)
lstmap = numpy.zeros(shape=(60*len(mslist)))
flagmap = numpy.zeros(shape=(6,1024,60*len(mslist)))
ms_lst_delta = 0.0937
ms_ts_delta = ms_lst_delta/60
ms = casac.casac.ms()
for i,msfile  in enumerate(mslist):


    
    print(i)
    ms.open(msfile)
    mssplit = msfile.split('.')
    ms_lst = mssplit[5]+'.'+mssplit[6]
    ms_lst = float(ms_lst)
    print(ms_lst)

    # Remove autocorrelations
    #ms.selecttaql("ANTENNA1!=ANTENNA2") 
    data = ms.getdata(['antenna1','antenna2','amplitude','phase','flag'],ifraxis=True)
    ms.done()
    amp = data['amplitude']
    phase = data['phase']
    ant1 = data['antenna1']
    ant2 = data['antenna2']
    flags = data['flag']
    ants = numpy.vstack((ant1,ant2))
    ants = ants.T
    dtype = 'int, int'    
    for l in numpy.arange(60):
        lstmap[i*60 + l] = ms_lst + l * ms_ts_delta
    ants = ants.astype(int).copy()
    ants = ants.view(dtype)
    print(numpy.shape(amp))
    for j, triads in enumerate(listoflistoftriads):
        print(j)
        triads = numpy.array(triads,dtype=int).copy()
        print("TRIAD NUMBER: ")
        print(numpy.shape(triads))
        triads = triads.view(dtype)
        #antv = ants.copy()

        fwd = numpy.in1d(ants,triads)
        print(fwd)
        bkwd = numpy.in1d(ants[:,::-1],triads)
        print(bkwd)
        #sign = fwd.astype(int) - bkwd.astype(int)
        value = fwd | bkwd
        
        ind = numpy.where(value)[0]
        
        print(ind.shape)
        #ind = ind * sign
        print(ind)

        amps = amp[0,:,ind,:]
        phases = phase[0,:,ind,:]
        print(numpy.shape(amps))

        comp = amps * numpy.exp(1j * phases)
        print("COMP SHAPE")
        print(comp.shape)
        datam = numpy.mean(comp,axis=0) # This should be a circular average.
        print(numpy.shape(datam))
        try:
            snrmap[j,:,i*60:(i+1)*60] = datam
            flagmap[j,:,i*60:(i+1)*60] = numpy.logical_or.reduce(flags[0,:,:,:],axis=1)
        except ValueError: #Last one isn't the same size..

            pass



numpy.savez("snrmap.npz",visamp=snrmap,last=lstmap,flag=flagmap)
		
