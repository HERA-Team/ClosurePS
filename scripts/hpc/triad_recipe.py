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


#processing specification
ProcessData=True
Pol="yy" #or "yy", etc
RefAnt="53" #as a string
#InTimes = [13298]
#InDay=2458098
inAntenna=[0,   1,   2,  11,  12,  13,  14,  23,  24,  25,  26,  27,  36,
        37,  38,  39,  40,  41,  50,  51,  52,  53,  54,  55,  65,  66,
        67,  68,  69,  70,  71,  82,  83,  84,  85,  87, 120, 121, 122,
       123, 124, 141, 142, 143]
inTriads=[[0,11,12],[0,1,12],[1,12,13],[1,2,13],[2,13,14],[11,23,24],[11,12,24],[12,24,25],[12,13,25],[13,25,26],[13,14,26],[14,26,27],[23,36,37],[23,24,37],[24,37,38],[24,25,38],[25,38,39],[25,26,39],[26,39,40],[26,27,40],[27,40,41],[36,37,51],[37,51,52],[37,38,52],[38,52,53],[38,39,53],[39,53,54],[39,40,54],[40,54,55],[40,41,55],[51,66,67],[51,52,67],[53,54,69],[54,69,70],[54,55,70],[55,70,71],[65,66,82],[66,82,83],[66,67,83],[67,83,84],[70,71,87],[121,122,141],[122,141,142],[122,123,142],[123,142,143],[123,124,143]]


IDRDays=[2458050,2458051,2458052,2458054,2458055,2458056,2458058,2458059]
NUMBER_OF_CORES = len(IDRDays)
    
#string builders for various filenames
#InData = glob.glob(os.path.join(DataDir,str(InDay),Pol,"*.uv")) #process the whole directory
#InData = [ os.path.join(DataDir, str(InDay), Pol,"zen." + str(InDay)+"." + str(x) + "." + str(Pol) + ".HH.uv") for x in InTimes]
InData=[]
for d in IDRDays:
    [InData.append(f) for f in glob.glob(os.path.join(DataDir,str(d),Pol,"*.uv"))]
InData = sorted(InData)
#get all uv files in the directory
#print InData
#sys.exit()

#InData=InData[1:2]
###############################


def makemsfile(fin,**kwargs):
    
    print(fin)
    hh = c.hf(makemsfile, fin)
    mm = repo.get(hh)
    print(mm)
    if mm:
        c.trc( "[Cached]", "makemsfile", fin)
    else:
        c.trc( "[Eval] ", "makemsfile", fin, " -> ", hh)
        UV = pyuvdata.UVData()
        UV.read_miriad(fin)
        UV.phase_to_time(Time(UV.time_array[0], format='jd', scale='utc'))
        tempf = repo.mktemp()
        os.remove(tempf)
        UV.write_uvfits(tempf, spoof_nonessential=True)
        if not os.path.exists(tempf):
            raise RuntimeError("No output produced by mkuvfits!")
        foms = c.importuvfits(tempf)
        os.remove(tempf)
        #flms = c.flagdata(foms,autocorr=True)
        mm = repo.put(foms, hh)
    return(mm)

def genclosurephase(fin,**kwargs):

    mm=makemsfile(fin,kwargs)
    fout=os.path.split(fin)[-1]+".npz"    
    r=hc.closurePh(mm,trlist=inTriads,alist=inAntenna)
    np.savez(fout,**r)
    if not os.path.exists(fout):
        raise RuntimeError("No output produced by hc.closurePh !")

    return(fout)

def genclosurephase(fin,**kwargs):
    print(fin)
    hh = c.hf(genclosurephase, fin)
    mm = repo.get(hh)
    print(mm)
    if mm:
        c.trc( "[Cached]", "makemsfile", fin)
    else:
        c.trc( "[Eval] ", "makemsfile", fin, " -> ", hh)
        UV = pyuvdata.UVData()
        UV.read_miriad(fin)
        UV.phase_to_time(Time(UV.time_array[0], format='jd', scale='utc'))
        tempf = repo.mktemp()
        os.remove(tempf)
        UV.write_uvfits(tempf, spoof_nonessential=True)
        if not os.path.exists(tempf):
            raise RuntimeError("No output produced by mkuvfits!")
        foms = c.importuvfits(tempf)
        os.remove(tempf)
        flms = c.flagdata(foms,autocorr=True)
        mm = repo.put(flms, hh)
    
    fout=os.path.split(fin)[-1]+".npz"    
    r=hc.closurePh(mm,trlist=inTriads,alist=inAntenna)
    np.savez(fout,**r)
    if not os.path.exists(fout):
        raise RuntimeError("No output produced by hc.closurePh !")
       



    return(fout)

	
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
