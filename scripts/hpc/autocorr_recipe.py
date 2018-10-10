import os,sys,glob
from recipe import casatasks as c
from recipe import repo
import casac
import casa
import shutil
import pyuvdata
import heracasa.closure as hc
import inspect
import numpy as np

#File location specification
DataDir="/rds/project/bn204/rds-bn204-asterics/HERA/data"
c.repo.REPODIR="/rds/project/bn204/rds-bn204-asterics/cache"

#processing specification
ProcessData=True
Pol="xx" #or "yy", etc
RefAnt="53" #as a string
#InTimes = [13298]
#InDay=2458098
inAntenna=[0,   1,   2,  11,  12,  13,  14,  23,  24,  25,  26,  27,  36,
        37,  38,  39,  40,  41,  50,  51,  52,  53,  54,  55,  65,  66,
        67,  68,  69,  70,  71,  82,  83,  84,  85,  87, 120, 121, 122,
       123, 124, 140, 141, 142, 143]
inTriads=[[0,11,12],
          [0,1,12],
          [1,12,13],
          [1,2,13],
          [2,13,14],
          [11,23,24],
          [11,12,24],
          [12,24,25],
          [12,13,25],
          [13,25,26],
          [13,14,26],
          [14,26,27],
          [23,36,37],
          [23,24,37],
          [24,37,38],
          [24,25,38],
          [25,38,39],
          [25,26,39],
          [26,39,40],
          [26,27,40],
          [27,40,41],
          [36,37,51],
          [37,51,52],
          [37,38,52],
          [38,52,53],
          [38,39,53],
          [39,53,54],
          [39,40,54],
          [40,54,55],
          [40,41,55],
          [51,66,67],
          [51,52,67],
          [53,54,69],
          [54,69,70],
          [54,55,70],
          [55,70,71],
          [65,66,82],
          [66,82,83],
          [66,67,83],
          [67,83,84],
          [121,140,141],
          [121,122,141],
          [122,141,142],
          [123,142,143],
          [123,124,143]]




# [[0,11,12],
#           [0,1,12],
#           [1,12,13],
#           [1,2,13],
#           [2,13,14],
#           [11,23,24],
#           [11,12,24],
#           [12,24,25],
#           [12,13,25],
#           [13,25,26],
#           [13,14,26],
#           [14,26,27],
#           [23,36,37],
#           [23,24,37],
#           [24,37,38],
#           [24,25,38],
#           [25,38,39],
#           [25,26,39],
#           [26,39,40],
#           [26,27,40],
#           [27,40,41],
#           [36,37,51],
#           [37,51,52],
#           [37,38,52],
#           [38,52,53],
#           [38,39,53],
#           [39,53,54],
#           [39,40,54],
#           [40,54,55],
#           [40,41,55],
#           [51,66,67],
#           [51,52,67],
#           [53,54,69],
#           [54,69,70],
#           [54,55,70],
#           [55,70,71],
#           [65,66,82],
#           [66,82,83],
#           [66,67,83],
#           [67,83,84],
#           [82,99,100],
#           [82,83,100],
#           [83,100,101],
#           [83,84,101],
#           [84,101,102],
#           [99,100,118],
#           [100,118,119],
#           [100,101,119],
#           [101,119,120],
#           [101,102,120],
#           [102,120,121],
#           [102,103,121],
#           [103,121,122],
#           [118,137,138],
#           [118,119,138],
#           [119,138,139],
#           [119,120,139],
#           [120,139,140],
#           [121,140,141],
#           [121,122,141],
#           [122,141,142],
#           [123,142,143],
#           [123,124,143]]

IDRDays=[2458098,2458099,2458101,2458102,2458103,2458104,2458105,2458106,2458107,2458108,2458109,2458110,2458111,2458113,2458114,2458115,2458116,2458118,2458119,2458120,2458121,2458122,2458123,2458124,2458125]
#string builders for various filenames
#InData = glob.glob(os.path.join(DataDir,str(InDay),Pol,"*.uv")) #process the whole directory
#InData = [ os.path.join(DataDir, str(InDay), Pol,"zen." + str(InDay)+"." + str(x) + "." + str(Pol) + ".HH.uv") for x in InTimes]
InData=[]
for d in IDRDays:
    [InData.append(f) for f in glob.glob(os.path.join(DataDir,str(d),Pol,"*.uv"))]
#get all uv files in the directory
#print InData
#sys.exit()

#InData=InData[1:2]
###############################
def mkuvfits(fin):
    hh=c.hf(mkuvfits, fin)
    mm=repo.get(hh)
    if mm:
        c.trc( "[Cached]", "mkuvfits", fin)
        return mm
    else:
        c.trc( "[Eval] ", "mkuvfits", fin, " -> ", hh)
        UV = pyuvdata.UVData()
        UV.read_miriad(fin,'miriad')
        UV.phase_to_time(UV.time_array[0])
        tempf=repo.mktemp()
        os.remove(tempf)
        UV.write_uvfits(tempf,'uvfits')
        if not os.path.exists(tempf):
            raise RuntimeError("No output produced by mkuvfits !")
        return repo.put(tempf, hh)

def extract_autocorrelations(fin, **kwargs):


    hh=c.hf(extract_autocorrelations, fin)
    mm=repo.get(hh)
    if mm:
        c.trc( "[Cached]", "mkuvfits", fin)
        return mm
    else:
        c.trc( "[Eval] ", "mkuvfits", fin, " -> ", hh)
        ms=casac.casac.ms()
        ms.open(fin)
        tempf=repo.mktemp()
        os.remove(tempf)
        data=ms.getdata(['antenna1','antenna2'.'data'],ifraxis=True)
        autocorrs_ind = numpy.arghwhere(data['antenna1'] == data['antenna2'])
        visrecords = data['data']
        autocorrs = visrecords[:,:,tuple(autocorrs_ind),:]
        np.savez(tempf,**autocorrs)
    if not os.path.exists(tempf+".npz"):
        raise RuntimeError("No output produced by autocorrelations :(")
    return repo.put(tempf+".npz",hh)
    
    

    
def mkclosurephase(fin, **kwargs):
    """Implement the closure phase calculation from heracasa with recipe"""
    hh=c.hf(mkclosurephase,inspect.getcallargs(hc.closurePh,fin,kwargs))
    mm=repo.get(hh)
    if mm:
        c.trc("[Cached]","hc.closurePh",fin,kwargs)
        return mm
    else:
        c.trc("[Eval]","hc.closurePh",fin,kwargs)
        tempf=repo.mktemp()
        os.remove(tempf)
        r=hc.closurePh(fin,**kwargs)
        np.savez(tempf,**r)
    if not os.path.exists(tempf+".npz"):
        raise RuntimeError("No output produced by hc.closurePh !")
    return repo.put(tempf+".npz",hh)

def copyfileoutput(dsin,dsout,extension):
    """copy output from dsin using dsout as directory names, clobber everything"""
    for i,f in enumerate(dsin):
        fin=f
        fout=os.path.split(dsout[i])[-1]+"."+extension
        if os.path.isfile(fout):
            os.remove(fout)
        shutil.copyfile(fin,fout)
    return fout
	
###############################
def main():
    f=[mkuvfits(fn) for fn in InData]
    fi = [c.importuvfits(fn) for fn in f]
    fif = [c.flagdata(fn, autocorr=True) for fn in fi]
    fifc = [extract_autocorrelations(fn) for fn in fif]
    print copyfileoutput(fifc,InData,"npz")
if ProcessData: main()

