# Whole of season 2 analysis

import numpy

if 0:
    fin=numpy.load("/home/bnikolic/data/HERA/closure/eq28_xx_binned.npz")

if 1:
    pyplot.clf()
    pyplot.plot(fin["last"][100], label="Record #100")
    pyplot.title("Time drift for aligned record across days" )
    pyplot.xlabel("Day number")
    pyplot.ylabel("Apparent LST")
    pyplot.savefig("plots/timedrift.png")
    
    
