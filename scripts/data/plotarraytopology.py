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
    
    current_rot_angle = numpy.arctan2(ant_y[-1],ant_x[-1])
    print("Current Rotation: ")
    print(current_rot_angle)
    
    # Rotate antenna positions
    angle_rot = -current_rot_angle
    cost = numpy.cos(angle_rot)
    sint = numpy.sin(angle_rot)
    
    xcost = ant_x * cost
    xsint = ant_x * sint
    ysint = ant_y * sint
    ycost = ant_y * cost
    
    x_rot = xcost - ysint
    y_rot = xsint + ycost
    

    fig, ax = plt.subplots(figsize=(8,8))
    ax.scatter(x_rot,y_rot)

    for i,label in enumerate(ant_nums):
        ax.annotate(label,(x_rot[i],y_rot[i]))

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.savefig("hera_array.pdf")
    plt.show()
    
if __name__=="__main__":
    main()
    
    
