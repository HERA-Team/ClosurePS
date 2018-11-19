import numpy
import matplotlib.pyplot as plt
import pyuvdata
import argparse
import textwrap



def main():

    command = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''
    ---------------------------------------------------------------------------------------
    HERA Array Topology Viewer

    Author: James Kent
    Institution: University of Cambridge, 2018
    Email: jck42@cam.ac.uk

    Takes a HERA miriad file and plots the array topology.

    Takes x,y locations to plot top-down view of the array, with antenna numbers.
'''))
    command.add_argument('filepath', help = ('miriad file from HERA'))
    args = command.parse_args()

    UV = pyuvdata.UVData()
    UV.read_miriad(args.filepath)
    ant_locs = UV.antenna_positions
    ant_nums = UV.antenna_numbers

    ant_x = ant_locs[:,0]
    ant_y = ant_locs[:,1]

    fig, ax = plt.subplots()
    ax.scatter(ant_x,ant_y)

    for i,label in enumerate(ant_nums):
        ax.annotate(label,(ant_x[i],ant_y[i]))

    plt.show()
    
if __name__=="__main__":
    main()
    
    
