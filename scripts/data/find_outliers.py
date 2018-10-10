import numpy, os
import scipy
import argparse
import matplotlib.pyplot as plt

def outlier_posterior(closure_chan,
                      concentration):

    prior = 1 / numpy.shape(closure_chan)[0]
    
    #for closure in closure_chan:



def plot_unit_circle(closures_chan):

    xc = numpy.cos(closures_chan)
    yc = numpy.sin(closures_chan)

    return xc, yc
    
def conc_measure(angular_data):
    
        
    x, y = plot_unit_circle(angular_data)
    concentration = numpy.sqrt(numpy.square(numpy.mean(x)) + numpy.square(numpy.mean(y)))
    return concentration

#def average_triads(closures):
#
#    for i in rang

## Takes closures in and finds the MSE for each triad with each other.
# Might allow grouping coherent antennas.

def calculate_pse_matrix(closures):

    pse_vals = numpy.ma.masked_array(numpy.zeros(shape=(numpy.shape(closures)[1],
                                                         numpy.shape(closures)[1])))

    closures[0,:,:150] = 0.0
    closures[0,:,800:] = 0.0
    cls_x = numpy.cos(closures)
    cls_y = numpy.sin(closures)
    ## More rowing bants: but with more cross-correlation
    for i in range(numpy.shape(closures)[1]):
        cldf_x = cls_x[0,i,:]
        cldf_y = cls_y[0,i,:]
        for j in range(numpy.shape(closures)[1]):
            pse_vals[i,j] = numpy.sum(numpy.sqrt(numpy.square(cldf_x - cls_x[0,j,:]) +
                                                 numpy.square(cldf_y - cls_y[0,j,:])))/(1024-150-(1024-800))
            
    print(pse_vals)
    fig = plt.figure(figsize=(6, 3.2))
    ax = fig.add_subplot(111)
    ax.set_title("Distance Metric Values")
    plt.imshow(pse_vals)
    plt.xlabel("Triad No.") 
    plt.ylabel("Triad No.")
    plt.colorbar(orientation='vertical')
    plt.show()

    
    for distance in pse_vals:

        sorted_arr = numpy.argsort(distance)
        print(sorted_arr)
        
    return pse_vals



## Takes triads in the [numpy.array] form of dimensionality [x,y,z] where:
# x = number of LST's to average over
# y = number of triads in saple
# z = number of channels


def plot_closures(closures):

    xvar = numpy.arange(0,1024) * (100.0/1024.0) + 100
    
    for i in range(0,numpy.shape(closures)[1]):
        #plt.plot(xvar, closures[0,i,:],'x')
        plt.plot(xvar, closures[0,i,:] + i * 2 * numpy.pi,',')
        plt.show(block=False)
        #x = input()
    plt.title("Closure Phases (offset by 2 * pi) at LST x.37538083174 on 2458098.58782")
    plt.ylim((-4,150))
    plt.ylabel("Phase (radians)")
    plt.xlabel("Frequency (MHz)")
    plt.show()

def brute_force_var(closures):

    std_devt = numpy.ma.masked_array(numpy.empty(shape=(numpy.shape(closures)[1],numpy.shape(closures)[2])))
    std_dev_dev = numpy.ma.masked_array(numpy.empty(shape=(numpy.shape(closures)[1], numpy.shape(closures)[2])))
    
    #First work out standard deviations with no triads masked out.
    cl = closures
    cl[:,:,819:] = 0.0
    cl[:,:,:150] = 0.0
    cx = numpy.cos(cl)
    cy = numpy.sin(cl)
    cxm = numpy.mean(cx, axis=1)
    cym = numpy.mean(cy, axis=1)

    conc_measure_all = numpy.empty(shape=numpy.shape(cxm)[1])
    for j in range(numpy.shape(cxm)[1]):
        conc_measure_all[j] = numpy.sqrt(numpy.square(numpy.mean(cxm[:,j]))
                                     + numpy.square(numpy.mean(cym[:,j])))
    std_dev_all = numpy.sqrt(-2 * numpy.log(conc_measure_all))
    
    ## Rowing banter: Seat race the triads, make them sweat to stay in the dataset.
    for i in range(numpy.shape(closures)[1]):
        cld = numpy.ma.masked_array(closures)
        cld[:,i,:] = numpy.ma.masked
        cld[:,:,819:] = 0.0
        cld[:,:,:150] = 0.0

        cld_x = numpy.cos(cld)
        cld_y = numpy.sin(cld)
        cld_x_m = numpy.ma.mean(cld_x, axis=1)
        cld_y_m = numpy.ma.mean(cld_y, axis=1)

        conc_measure = numpy.empty(shape=numpy.shape(cld_x_m)[1])
        for j in range(numpy.shape(cld_x_m)[1]):
            conc_measure[j] = numpy.sqrt(numpy.square(numpy.mean(cld_x_m[:,j])) + numpy.square(numpy.mean(cld_y_m[:,j])))

        std_dev = numpy.sqrt(-2 * numpy.log(conc_measure))
        #print(numpy.shape(conc_measure))
        #print(numpy.shape(std_dev))
        #print(numpy.shape(std_devt[i,:]))
        std_dev[numpy.isnan(std_dev)] = 0.0
        std_devt[i,:] = std_dev
        std_dev_dev[i,:] = std_dev - std_dev_all
        #print(std_dev[510])
        #print(numpy.shape(cld_x_m))


    ## Find sections where the standard deviation decreases as a result of removing that triad. 
    std_dev_find_bad_bins = numpy.zeros(shape=(numpy.shape(closures)[1],numpy.shape(closures)[2]))
    for i in range(numpy.shape(closures)[2]):

        std_dev_find_bad_bins[:,i] = (std_devt[:,i] < std_dev_all[i])

    sum_of_chan = numpy.sum(std_dev_find_bad_bins,axis=1)
    print(sum_of_chan)

    xvar = numpy.arange(0,1024) * (100.0/1024.0) + 100

    triads = range(numpy.shape(closures)[1])
    for i in triads:
        plt.plot(xvar, std_dev_dev[i,:],'x')
        plt.show(block=False)
        print(i)
        x = raw_input()

    plt.title("Triad Seat Racing: Deviation from Mean Std. Dev")
    plt.ylabel("Deviation from Mean")
    plt.xlabel("Frequency(MHz)")
    plt.show()

    #for i in triads:
    #    plt.plot(xvar,closures(

        
        
    # print("Average Std Devs")
    # xvar = numpy.arange(0,1024)
    # for i in range(numpy.shape(closures)[1]):
    #     print(numpy.ma.mean(std_devt[i,:]))
    #     plt.plot(xvar, std_devt[i,:],'x')
    #     plt.show(block=False)
    #     x = input()
    # plt.show()
    
def main():

    args = None
    command = argparse.ArgumentParser()
    command.add_argument('filepath', help = ('.npz directory'))
    command.add_argument('working_directory', help = ('Working Directory'))
    #command.add_argument('flag_file',help=('Flag file for channels'))
    args = command.parse_args()

    if (os.path.exists(args.working_directory) == False):
        os.makedirs(args.working_directory)
    
    files = os.listdir(args.filepath)
    npz_files = []
    for file in files:
        if file.endswith(".npz"):
                    npz_files.append(file)

    npz_array = numpy.empty(shape=(len(npz_files),26,1024))
    for i,npz in enumerate(npz_files):
        with numpy.load(args.filepath+npz) as data:
            npz_array[i,:,:] = data.f.arr_0
#            numpy.append(npz_array,data.f.arr_0)

    print(numpy.shape(npz_array))
    #print(npz_array[4,13,:])
    plot_closures(npz_array)
    #brute_force_var(npz_array)
    calculate_pse_matrix(npz_array)

if __name__=="__main__":
    main()
