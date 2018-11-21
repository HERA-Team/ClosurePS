
import os
import numpy
import argparse
import matplotlib.pyplot as plt

def main():


    # Command Line Parse:
    
    args = None
    command = argparse.ArgumentParser()
    command.add_argument('filepath', help = ('.npz directory'))
    command.add_argument('working_directory', help = ('Working Directory. Where files will be put.'))
    command.add_argument('--multiple_days', action='store_true',required=False, help=('Will parse a directory of closures from different days and sort through them'))
    args = command.parse_args()

    if (os.path.exists(args.working_directory) == False):
        os.makedirs(args.working_directory)

    files = os.listdir(args.filepath)
    npz_files =[]
    for file in files:
        if file.endswith(".npz"):
            npz_files.append(args.filepath+file)


    if args.multiple_days:
        detected_datestamps = []
        for file in npz_files:
            detected_datestamps.append(file.split('.')[1])
        detected_datestamps = set(detected_datestamps)
        detected_datestamps = sorted(detected_datestamps)

        
        #First find all files with a certain datestamp
        for datestamp in detected_datestamps:
            all_closures = None
            all_last = None
            for file in npz_files:
                if file.split('.')[1] == datestamp:
                    with numpy.load(file) as data:
                        closures = data['phase']
                        last=data['LAST']
                        all_closures = numpy.zeros(shape=numpy.shape(closures))
                        all_last = numpy.zeros(shape=numpy.shape(last))
                    break
            for file in sorted(npz_files):
                if file.split('.')[1] == datestamp:
                    with numpy.load(file) as data:
                        closures = data['phase']
                        flags = data['flags']
                        last = data['LAST']
                        all_closures=numpy.concatenate((all_closures,closures),axis=3)
                        all_last=numpy.concatenate((all_last,last))
            all_closures = all_closures[:,:,:,60:]
            all_last = all_last[60:]
            all_closures_x = numpy.cos(all_closures)
            all_closures_y = numpy.sin(all_closures)
            averaged_x = numpy.mean(all_closures_x, axis=0)
            averaged_y = numpy.mean(all_closures_y, axis=0)
            r = numpy.sqrt(numpy.square(averaged_x) + numpy.square(averaged_y))
            std_dev = numpy.sqrt(-2 * numpy.log(r))

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title("LST Quality on JD "+datestamp)



            plt.imshow(std_dev[0,:,:])


            all_lsth = [(lst%1) * 24 for lst in all_last]
            all_lsth = numpy.asarray(all_lsth)
            integer_ticks = []
            for i,tick in enumerate(all_lsth):
                if tick - numpy.floor(tick) < 0.003:
                    integer_ticks.append(i)
                    all_lsth[i] = numpy.floor(tick)
            xticks = ["{:2.1f}".format(lst) for lst in all_lsth]
            ax.set_xticks(numpy.arange(0,numpy.shape(all_last)[0]))
            ax.set_xticklabels(xticks)
            for label in ax.xaxis.get_ticklabels():
                label.set_visible(False)                
            for i in integer_ticks:
                ax.axvline(i, linestyle='--', color='k')
                label = ax.xaxis.get_ticklabels()[i]
                label.set_visible(True)
            
                    #ax.set_xticklabels(["{:6.4f}".format(i) for i in all_last])
            plt.ylabel("Channel")
            plt.xlabel("Local Sidereal Time")
            # ax.set_xlabel((all_last[0],all_last[-1]))
            plt.clim(0,2.5)
            plt.colorbar(orientation='horizontal')
            plt.savefig(args.working_directory+datestamp+".eps")

            

    # Case of a directory with just one file in it, or concatenate EVERY day together. Why would you do this?
    else:
            
        # Concatenate all the .npz files:
        all_closures = None
        all_last = None
        for file in npz_files:
            with numpy.load(file) as data:
                closures = data['phase']
                last = data['LAST']
                all_closures = numpy.zeros(shape=numpy.shape(closures))
                all_last = numpy.zeros(shape=numpy.shape(last))
            break
        print(numpy.shape(all_closures))
        print(numpy.shape(all_last))
        for file in sorted(npz_files):
            print(file)
            with numpy.load(file) as data:
                closures = data['phase']
                flags = data['flags']
                last = data['LAST']
                print(numpy.shape(closures))
                all_closures = numpy.concatenate((all_closures,closures),axis=3)
                all_last = numpy.concatenate((all_last, last))

        # Find standard deviation across triads at every timestamp and channel:

        all_closures = all_closures[:,:,:,60:]
        all_last = all_last[60:]
        all_closures_x = numpy.cos(all_closures)
        all_closures_y = numpy.sin(all_closures)
    
        averaged_x = numpy.mean(all_closures_x,axis=0)
        averaged_y = numpy.mean(all_closures_y,axis=0)
        r = numpy.sqrt(numpy.square(averaged_x) + numpy.square(averaged_y))
        std_dev = numpy.sqrt(-2 * numpy.log(r))


    
        print(numpy.argwhere(numpy.logical_and(all_last < 64834.1067, all_last > 64834.08)))
        print(numpy.shape(numpy.logical_and(all_last < 64834.1067, all_last > 64834.08)))
    
        print("All Closures/Std_dev Shape:")
        print(numpy.shape(all_closures))
        print(numpy.shape(std_dev))
        print("----------------------------")

    

    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title("LST Quality on JD 2458160")



        plt.imshow(std_dev[0,:,:])


        all_lsth = [(lst%1) * 24 for lst in all_last]
        all_lsth = numpy.asarray(all_lsth)
        integer_ticks = []
        for i,tick in enumerate(all_lsth):
            if tick - numpy.floor(tick) < 0.003:
                integer_ticks.append(i)
                all_lsth[i] = numpy.floor(tick)
        xticks = ["{:2.1f}".format(lst) for lst in all_lsth]
        ax.set_xticks(numpy.arange(0,numpy.shape(all_last)[0]))
        ax.set_xticklabels(xticks)
        for label in ax.xaxis.get_ticklabels():
            label.set_visible(False)                
        for i in integer_ticks:
            ax.axvline(i, linestyle='--', color='k')
            label = ax.xaxis.get_ticklabels()[i]
            label.set_visible(True)
           
        plt.ylabel("Channel")
        plt.xlabel("Local Sidereal Time")
        # ax.set_xlabel((all_last[0],all_last[-1]))
        plt.clim(0,2.5)
        plt.colorbar(orientation='horizontal')
        plt.show()

        
if __name__=="__main__":
    main()
