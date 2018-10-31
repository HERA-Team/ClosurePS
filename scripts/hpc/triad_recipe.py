import os,sys,glob
from recipe import casatasks as c
from recipe import repo
import casac
import casa
import shutil
import pyuvdata
import heracasa.closure as hc
import inspect
import multiprocessing
import numpy as np

from astropy.time import Time
#File location specification
DataDir="/rds/project/bn204/rds-bn204-asterics/HERA/data"
c.repo.REPODIR="/rds/project/bn204/rds-bn204-asterics/cache"

#Multiprocessing specification
NUMBER_OF_CORES = 16

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
inTriads=[[0,11,12],[0,1,12],[1,12,13],[1,2,13],[2,13,14],[11,23,24],[11,12,24],[12,24,25],[12,13,25],[13,25,26],[13,14,26],[14,26,27],[23,36,37],[23,24,37],[24,37,38],[24,25,38],[25,38,39],[25,26,39],[26,39,40],[26,27,40],[27,40,41],[36,37,51],[37,51,52],[37,38,52],[38,52,53],[38,39,53],[39,53,54],[39,40,54],[40,54,55],[40,41,55],[51,66,67],[51,52,67],[53,54,69],[54,69,70],[54,55,70],[55,70,71],[65,66,82],[66,82,83],[66,67,83],[67,83,84],[70,71,87],[120,121,140],[121,140,141],[121,122,141],[122,141,142],[122,123,142],[123,142,143],[123,124,143]]

IDRDays=[2458141,2458142,2458143,2458144,2458145,2458146,2458147,2458148,2458149,2458150,2458151,2458152,2458153,2458154,2458155,2458156,2458157,2458158,2458159,2458160]

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


def genclosurephase(fin,**kwargs):
    hh = c.hf(mkuvfits, fin)
    mm = repo.get(hh)
    print(mm)
    if mm:
        c.trc( "[Cached]", "mkuvfits", fin)
    else:
        c.trc( "[Eval] ", "mkuvfits", fin, " -> ", hh)
        UV = pyuvdata.UVData()
        UV.read_miriad(fin)
        UV.phase_to_time(Time(UV.time_array[0], format='jd', scale='utc'))
        tempf = repo.mktemp()
        os.remove(tempf)
        UV.write_uvfits(tempf, spoof_nonessential=True)
        if not os.path.exists(tempf):
            raise RuntimeError("No output produced by mkuvfits!")
        
        mm = repo.put(tempf, hh)

    foms = c.importuvfits(mm)
    flms = c.flagdata(foms,autocorr=True)
    
    hh = c.hf(mkclosurephase, inspect.getcallargs(hc.closurePh,flms,trlist=inTriads,alist=inAntenna))
    mm=repo.get(hh)
    print(mm)
    if mm:
        c.trc("[Cached]", "hc.closurePh",flms,kwargs)
    else:
        c.trc("[Eval]", "hc.closurePh",flms,kwargs)
        tempf=repo.mktemp()
        os.remove(tempf)
        r=hc.closurePh(flms,trlist=inTriads,alist=inAntenna)
        np.savez(tempf,**r)
        if not os.path.exists(tempf+".npz"):
            raise RuntimeError("No output produced by hc.closurePh !")
        mm = repo.put(tempf+".npz",hh)

    fout=os.path.split(fin)[-1]+".npz"
    if os.path.isfile(fout):
        os.remove(fout)
    shutil.copyfile(mm,fout)
    return(fout)

def mkuvfits(fin):
    hh=c.hf(mkuvfits, fin)
    mm=repo.get(hh)
    if mm:
        c.trc( "[Cached]", "mkuvfits", fin)
        return mm
    else:
        c.trc( "[Eval] ", "mkuvfits", fin, " -> ", hh)
        print(fin)
        UV = pyuvdata.UVData()
        UV.read_miriad(fin)
        UV.phase_to_time(Time(UV.time_array[0], format='jd', scale='utc'))
        tempf=repo.mktemp()
        os.remove(tempf)
        UV.write_uvfits(tempf,spoof_nonessential=True)
        if not os.path.exists(tempf):
            raise RuntimeError("No output produced by mkuvfits !")
        return repo.put(tempf, hh)

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


def multiprocess_wrapper(files,d):

    for f in files:
        genclosurephase(f,trlist=inTriads,alist=inAntenna)



def main():


    
    files_per_core = len(InData) / NUMBER_OF_CORES
    split_files = [InData[i:i+files_per_core] for i in range(0, len(InData), files_per_core)]
    remainder_files = InData[files_per_core * NUMBER_OF_CORES:len(InData)]

    jobs = []

    for i, list_slice in enumerate(split_files):
        print(list_slice)
        j = multiprocessing.Process(target=multiprocess_wrapper, args=(list_slice,'a'))
        jobs.append(j)

    for job in jobs:
        job.start()

    for file in remainder_files:
        genclosurephase(file,trlist=inTriads,alist=inAntenna)
    
#    f=[genclosurephase(fn,trlist=inTriads,alist=inAntenna) for fn in InData]
#    fi = [c.importuvfits(fn) for fn in f]
#    fif = [c.flagdata(fn, autocorr=True) for fn in fi]
#    fifc = [mkclosurephase(fn,trlist=inTriads,alist=inAntenna) for fn in fif]
    print f
if ProcessData: main()
