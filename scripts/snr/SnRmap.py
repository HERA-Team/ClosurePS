import numpy
import os
import glob
import casac

n_channels = 1024

mslist = glob.glob("*ms")
mslist = sorted(mslist)
print(mslist)

snrmap = numpy.zeros(shape=(1024,60*len(mslist)))
stddevmap = numpy.zeros(shape=(1024,60*len(mslist)))

ms = casac.casac.ms()
for i,msfile  in enumerate(mslist):
	print(i)
	ms.open(msfile)
	# Remove autocorrelations
	ms.selecttaql("ANTENNA1!=ANTENNA2") 
	data = ms.getdata(['amplitude'],ifraxis=True)
	ms.done()
	data = data['amplitude']
	datam = numpy.mean(data,axis=2)
	datas = numpy.std(data,axis=2)
	try:
		snrmap[:,i*60:(i+1)*60] = datam[0,:,:]
		stddevmap[:,i*60:(i+1)*60] = datas[0,:,:]
	except ValueError: #Last one isn't the same size..
		pass



numpy.savez("snrmap.npz",visamp=snrmap,stddev=stddevmap)
		
