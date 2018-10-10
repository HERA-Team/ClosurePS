import numpy
import os
import glob
import pyuvdata
import casac

import multiprocessing

Datadir = "/lustre/aoc/projects/hera/H1C_IDR2/IDR2_1/LSTBIN/one_group/grp1"
Outputdir = "/lustre/aoc/projects/hera/jkent/HERA_SnR/"
pol = "xx"

def mkuvfits(file,outputdir):
	
	head,tail = os.path.split(file)
	
	if os.path.exists(file+"fits") == False:
		UV =pyuvdata.UVData()
		UV.read_miriad(file,'miriad')
		UV.phase_to_time(UV.time_array[0])
		UV.write_uvfits(outputdir+tail+"fits",'uvfits')
	return file+"fits"

def multiprocess_wrapper(splitfiles):
	
	[mkuvfits(f,Outputdir) for f in splitfiles]

def main():

	# Import everything.
	listuv = glob.glob(os.path.join(Datadir,"*xx.LST*.uvOCRSL"))
	listuv = sorted(listuv)
	print(listuv)
	UVFits = []
	#[UVFits.append(mkuvfits(f,Outputdir)) for f in listuv]
	#[casa.importuvfits(f,f+".ms") for f in UVFits]	
	
	numberofcores = 8

	filespercore = len(listuv)/numberofcores
	splitfiles = [listuv[i:i+filespercore] for i in range(0, len(listuv),filespercore)]
	remainderfiles = listuv[filespercore* numberofcores: len(listuv)]

	jobs = []

	for i,filelist in enumerate(splitfiles):
		print(filelist)
		j = multiprocessing.Process(target=multiprocess_wrapper, args=(filelist,))
		jobs.append(j)

	for job in jobs:
		job.start()
	for file in remainderfiles:
		mkuvfits(file,OutputDir)	

if __name__ == "__main__":
	main()
