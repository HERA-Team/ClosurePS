### HERA LST Aligner
#
# James Kent, University of Cambridge
# jck42@cam.ac.uk
#
# This script takes a directory of closure files outputted from recipe using
# heracasa (Bojan Nikolic, bn204@cam.ac.uk), and from there bins them so that
# closure measurements over successive days are aligned in LST.

from __future__ import print_function

from astropy.time import Time
import matplotlib.pyplot as plt
import os
import sys
import numpy
import argparse
import textwrap

LST_DP = 6


class LST_Alignment(object):
    """
    Class takes a dictionary of datestamps and aligns them to each other.
    We take advantage of sidereal time moving 235.9 seconds per day with respect
    to the solar day.

    Attributes:

    closure_directory      [String] Base directory where heracasa .npz files are located.
    
    date_set               [List] Ordered dates , earliest -> latest
    
    date_dict              [Dictionary] Dictionary of all dates with filenames.

    timestamp_delta        [Float] Delta between sidereal day and solar day.

    delta_ind              [Int] How many timestamp indices drift per sidereal day.

    delta_rem              [Float] Remainder from delta_ind calculation.

    Member Functions:

    __init__()             Initialises an instance of class LST_Alignment.

    __extract_closures()   Opens all the .npz files and concatenates them into a single dictionary 
                           keyed by date/lst.

    __align_closures()     Aligns the closures by LST.

    align_timestamps()     Public function which builds the date_dict and aligns the LST's, and 
                           returns them.

    TODO: Implement destructive alignment.
    """

    def __init__(self,
                 closure_directory, 
                 ordered_dates, 
                 date_dict,
                 lst_start,
                 lst_end,
                 integration_time=10.7374, 
                 sidereal_delta=235.909,
                 destructive = True): 
        """
        Initialises the LST_Alignment Class which aligns the LST's over successive
        Julian Dates"

        Inputs:

        closure_directory   [String] Root directory of heracasa .npz files.
        
        ordered_dates       [List] All datestamps in order. Earliest -> Latest

        date_dict           [Dictionary] of all dates with filenames.

        integration_time    [Float] HERA Integration Time. Default = 10.7s.

        sidereal_delta      [Float] Delta between sidereal day and solar day.
        
        destructive         [Bool] Destructive alignment or not? [NOT IMPLEMENTED (YET)]
        """
        self.closure_directory = closure_directory
        self.date_set = ordered_dates
        self.date_dict = date_dict
        self.integration_time = integration_time
        self.timestamp_delta = sidereal_delta / self.integration_time
        self.delta_ind = round(self.timestamp_delta) #Get closest indice
        self.delta_rem = self.timestamp_delta % 1
        self.triad_no = 0

        self.lst_start = lst_start
        self.lst_end = lst_end

        # We need to tweak lst_end to make sure we get the same number
        # of LST records each day.

        self.remainder = (((self.lst_end - self.lst_start)*3600.)/self.integration_time)%1
        self.lst_end = self.lst_end - (self.remainder * self.integration_time)/3600.


    def __extract_closures(self):

        """
        Opens all of the .npz files in date_dict, and sorts them into a new dictionary
        which is keyed by both date and LST. Thus the key mechanism is as so:

        |
        | - Date_1 - 'phase' - [Numpy Array]
        |            'flags' - [Numpy Array]
        |            'triads' - [Numpy Array]
        | - Date_2 - 'phase' - [Numpy Array] 
        |            'flags' - [Numpy Array] 
        |            'triads' - [Numpy Array]
        | - etc.

        """

        closure_dict = {}
     #Only works in Python 2.x
        for date, npz_files in sorted(self.date_dict.iteritems()):
            print(".",end="")
            sys.stdout.flush()
            closure_dict[date] = {}

            last_concat = None
            phases_concat = None
            flags_concat = None
            triads_ar = None
            #print(npz_files)
            # This is based on the that npz_files is sorted CORRECTLY
            # in chronological order
            for npz_file in npz_files:
                with numpy.load(self.closure_directory + npz_file) as data:

                    # LST's is used to build the keys of the second tree layer,
                    # as we want to organise by LST for our alignment.
                    last = data['LAST']                    
                    last = ((last%1)*24).round(decimals=LST_DP)
                    phases = data['phase']
                    if self.triad_no == 0:
                        self.triad_no = phases.shape[0]
                    flags = data['flags']
                    triads_ar = data['tr']

                    try:
                        last_concat = numpy.concatenate((last_concat,last),axis=0)
                        phases_concat = numpy.concatenate((phases_concat,phases),axis=3)
                        flags_concat = numpy.concatenate((flags_concat,flags),axis=3)
                    except ValueError:
                        last_concat = last
                        phases_concat = phases
                        flags_concat = flags
            
            #print(phases_concat.shape)
            #print(flags_concat.shape)
            #print(last_concat.shape)
            closure_dict[date]['phase'] = phases_concat
            closure_dict[date]['flags'] = flags_concat
            closure_dict[date]['last'] = last_concat
            closure_dict[date]['triads'] = triads_ar
        #sys.exit(1)           
        return(closure_dict)

    def __align_closures(self, reference_lst, closures):

        """
        Does the alignment of the LST's over successive Julian Days. 
        Little bit shakey and basic but works okay for now.

        Inputs:

        reference_lst  [Dictionary] First LST in dataset for reference. 

        closures       [Dictionary] Dictionary of all closures.


        TODO: Tidy this up, reference_lst not needed?
        """
        print("Generating LST Array.. ",end='')
        #Generate Numpy array from closure_dictionary. Makes life significantly easier.

        initial_date = sorted(closures.iteritems())[0]

        reference_lst_array = reference_lst['last']
        reference_phase_array = reference_lst['phase']
        
#        ib = numpy.searchsorted(reference_lst_array,self.lst_start)
#        ie = numpy.searchsorted(reference_lst_array,self.lst_end)

        #We are searching a rotated array due to LST wrapping around at 24 hours.
        zeroth_lst_index = numpy.argmin(reference_lst_array)        
        #Make sure it's 0 hrs LST
        if(reference_lst_array[zeroth_lst_index] < 0.0035):
            last_lower = reference_lst_array[:zeroth_lst_index]
            last_upper = reference_lst_array[zeroth_lst_index:]
                
            # Where we are aligning over zero lst, has to be carefully dealt with
            # by using a standard divide and conquer approach
            if(self.lst_start > self.lst_end):
                ib = numpy.searchsorted(last_lower,self.lst_start)
                ie = numpy.searchsorted(last_upper,self.lst_end)+zeroth_lst_index
                
#                if (ib == reference_lst_array.shape[0]) or (ib == 0):
#                        print(last_lower)

                    
            else:
                ib = numpy.searchsorted(last_lower,self.lst_start)
                
                if ib == 0:
                    ib = numpy.searchsorted(last_upper,self.lst_start)+zeroth_lst_index
#                    print(ib,end=" ")
                    ie = numpy.searchsorted(last_upper,self.lst_end)+zeroth_lst_index
#                    print(ie)

        else:
            ib = numpy.searchsorted(reference_lst_array,self.lst_start)
#            print(ib,end=" ")
            ie = numpy.searchsorted(reference_lst_array,self.lst_end)
#            print(ie)
        
        last_length = ie-ib
        dates = sum([1 for date in closures.iteritems()])
        triad_no = reference_phase_array.shape[0]
        channel_no = reference_phase_array.shape[2]
        days_arr = numpy.asarray([date[0] for date in sorted(closures.iteritems())])
        phases_subset_arr = numpy.zeros(shape=(triad_no,dates,channel_no,last_length))
        flags_subset_arr = numpy.zeros(shape=(triad_no,dates,channel_no,last_length))
        triads_arr = reference_lst['triads']
        lst_subset_arr = numpy.zeros(shape=(dates,last_length))
        bad_days = []

        
        for i,date in enumerate(sorted(closures.iteritems())):
            print(".",end="")
            #print(date)
            datel = date[0]
            last = closures[datel]['last']

            is_sorted = numpy.logical_or.reduce(numpy.diff(last)>=0)
#            print(is_sorted)            
            phases = closures[datel]['phase']
            flags = closures[datel]['flags']

            #We are searching a rotated array due to LST wrapping around at 24 hours.
            zeroth_lst_index = numpy.argmin(last)            
            #Make sure its 0 hrs LST
            if(last[zeroth_lst_index] < 0.0035):
                last_lower = last[:zeroth_lst_index]
                last_upper = last[zeroth_lst_index:]
                
                # Where we are aligning over zero lst, has to be carefully dealt with
                # by using a standard divide and conquer approach
                if(self.lst_start > self.lst_end):
                    date_lst_lb = numpy.searchsorted(last_lower,self.lst_start)
                    date_lst_ub = numpy.searchsorted(last_upper,self.lst_end)+zeroth_lst_index

                    if (date_lst_lb == last.shape[0]) or (date_lst_lb == 0):
                        print(last_lower)

                    
                else:
                    date_lst_lb = numpy.searchsorted(last_lower,self.lst_start)

                    if date_lst_lb == 0:
                        date_lst_lb = numpy.searchsorted(last_upper,self.lst_start)+zeroth_lst_index
                        date_lst_ub = numpy.searchsorted(last_upper,self.lst_end)+zeroth_lst_index

            else:
                date_lst_lb = numpy.searchsorted(last,self.lst_start)
                date_lst_ub = numpy.searchsorted(last,self.lst_end)
            # This looks messy but basically we just see if a day is within one record
            # if not we say its a bad day and get rid of it.
            try:
                phases_subset_arr[:,i,:,:] = phases[:,0,:,date_lst_lb:date_lst_ub]
                flags_subset_arr[:,i,:,:] = flags[:,0,:,date_lst_lb:date_lst_ub]
                lst_subset_arr[i,:] = last[date_lst_lb:date_lst_ub]

            except ValueError:
                try:
                    phases_subset_arr[:,i,:,:] = phases[:,0,:,date_lst_lb:date_lst_ub-1]
                    flags_subset_arr[:,i,:,:] = flags[:,0,:,date_lst_lb:date_lst_ub-1]
                    lst_subset_arr[i,:] = last[date_lst_lb:date_lst_ub-1]
                except ValueError:

                    try:
                        phases_subset_arr[:,i,:,:] = phases[:,0,:,date_lst_lb:date_lst_ub+1]
                        flags_subset_arr[:,i,:,:] = flags[:,0,:,date_lst_lb:date_lst_ub+1]
                        lst_subset_arr[i,:] = last[date_lst_lb:date_lst_ub+1]
                    except ValueError: # At this point we just assume the day's data is corrupted or unlike the others.
#                        print("Bad day detected: ",datel)
                        bad_days.append(i)
                        continue

        if len(bad_days) > 0:
            bad_days = numpy.asarray(bad_days)
            bad_days_list = [sorted(closures.iteritems())[bad_day][0] for bad_day in bad_days]
            phases_subset_arr = numpy.delete(phases_subset_arr,bad_days,axis=1)
            flags_subset_arr = numpy.delete(flags_subset_arr,bad_days,axis=1)
            lst_subset_arr = numpy.delete(lst_subset_arr,bad_days,axis=0)
            days_arr = numpy.delete(days_arr,bad_days,axis=0)
        print("done")
        if len(bad_days) > 0:
            print("Bad Days detected: ",bad_days_list)
            print("Final Days: ")
            print(days_arr)
        
        self.flags = numpy.transpose(flags_subset_arr,(3,1,0,2))
        self.phases = numpy.transpose(phases_subset_arr,(3,1,0,2))
        self.last = lst_subset_arr
        self.triads = triads_arr
        self.days = days_arr

        print("Final array shapes: ")
        print(self.last.shape)
        print(self.phases.shape)
        print(self.flags.shape)
        
        
    def align_timestamps(self):
        """
        Generates the closure_dictionary and generates a numpy array of aligned
        LST's, which can be used by the LST_Binner class to extract closures of
        choice across successive days.
        """
        print("Aggregating closure dictionary to array (can take a while)...",end="")
        closure_dict = self.__extract_closures()
        print("done")
        #print closure_dict

#        print("LST Reference: ")
        lst_ref = self.date_set[0]
#        print(lst_ref)
        lst_ref = closure_dict[lst_ref] #We interpolate to the LST's from the first day of observations.
        #print(lst_ref)
        print(len(lst_ref))
        print("Performing alignment...", end="")
        self.__align_closures(lst_ref, closure_dict)
        print("done")
        

    def save_binned_lsts(self,filename):
        
        """
        Takes our binned closures and outputs them to a .npz file.

        Inputs:

        filename    [String] Filename to save .npz file.
        """

        if self.days is not None:
            numpy.savez(filename,
                        days=self.days,
                        last=self.last,
                        closures=self.phases,
                        flags=self.flags,
                        triads=self.triads)
        else:
            raise ValueError("LST's not binned")


#This parses all of the files and breaks them up into datestamps. Each datestamp is then aligned correctly. 
class Julian_Parse(object):
    """ 
    Class to take a set of filepaths spat out from heracasa and create a dictionary of all the
    files, keyed by their Julian Date.

    Also returns the set of all dates found in the directory.

    Attributes:

    filepaths              [List] All filepaths from the directory ending in .npz. Assumed to be heracasa.
    
    file_dict              [Dictionary] Dictionary of npz files, keyed by their date.

    date_set               [Set] The set of dates.
    
    Member Functions:

    __init__()             Initialises an instance of class Julian_Parse.

    __find_unique_dates()  Finds all unique dates and returns the ordered list of the set of dates.

    __build_closure_dict() Takes the ordered list of the set of dates, and filepaths and sorts 
                           them into file_dict.

    break_up_datestamps()  Creates file_dict.

    return_datestamps()    Returns self.file_dict, self.date_set.

    TODO: Some sort of date ranging so you can control how much data we push through the binner.
    """


    
    # We assume all .npz files have 60 timestamps in them.
    def __init__(self, filepaths, date_start, date_end, exclude_days,npz_size = 60):
 
        """
        Initialises the Julian Parse Class which takes a full directory of 
        heracasa .npz files and splits them by date.

        """

        self.filepaths = filepaths
        self.date_start = date_start
        self.date_end = date_end
        self.exclude_days = exclude_days
        self.file_dict = None
        self.date_set = None

    # Parses all datestamps and converts to a set (finds unique values)
    def __find_unique_dates(self, filepaths):

        """
        Find unique dates from a directory of heracasa files.

        Inputs:

        filepaths [List] List of all .npz filepaths

        """

        detected_datestamps = []

        for file in self.filepaths:
            detected_datestamps.append(file.split('.')[1])
        detected_datestamps = set(detected_datestamps) # Convert list to set.
        detected_datestamps = sorted(detected_datestamps) # Convert set to ordered list.
        return(detected_datestamps)

    # Takes set of datestamps and filepaths and sorts them into a dict.
    # Dict layout is [Key: datestamp, Data: list of all files in that datestamp]
    def __build_closure_dict(self, datestamps, filepaths):

        """
        Splits the filepaths into a dictionary keyed by their Julian Date

        Inputs:

        datestamps [List] An ordered list of the datestamps.

        filepaths [List] List of al .npz filepaths.

        """
        file_dict = {}

        for datestamp in datestamps:
            file_dict[datestamp] = [] #Empty list.

        for datestamp in datestamps:

            for file in filepaths:
                if file.split('.')[1] == datestamp:
                    file_dict[datestamp].append(file)
            file_dict[datestamp] = sorted(file_dict[datestamp]) #Holy shit this just works? How awesome is that.
        #print(file_dict)
        return(file_dict)


    
    # This builds a dictionary of all unique datestamps and their .npz files
    def break_up_datestamps(self):

        """
        Breaks up the directory of .npz files into a dictionary, keyed by date.

        """
        print("Parsing closure directory... ", end="")       
        detected_datestamps = self.__find_unique_dates(self.filepaths)
        print("done")
        print("Discovering dates within specified range...", end="")
        detected_datestamps = sorted(list(filter(lambda el: int(el) >= self.date_start and int(el) <= self.date_end, detected_datestamps)))
        for excluded_date in self.exclude_days:
            print(excluded_date)
            detected_datestamps = filter(lambda el: int(el) != int(excluded_date), detected_datestamps)
        print("done")
        print("Detected Dates: ")
        print(detected_datestamps)
        print("Building dictionary of dates and filepaths...", end="")
        file_dict = self.__build_closure_dict(detected_datestamps, self.filepaths)
        print("done")
        self.file_dict = file_dict
        self.date_set = detected_datestamps

    # Returns dictionary.
    def return_datestamps(self):

        """
        Returns file_dict and date_set

        """
        return self.file_dict, self.date_set
    


def main():


    # Parse directory for input/output.
    command = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''
    ---------------------------------------------------------------------------------------
    HERA LST Aligner and Binner

    Author: James Kent
    Institution: University of Cambridge, 2018
    Email: jck42@cam.ac.uk

    Takes a directory of heracasa created closure phases (best created using recipe), 
    and concatenates all of the .npz files together, aligns the sidereal times, and 
    outputs and LST range of your choice.

    Can also optinally average closures across LSTs as well as calculate the standard
    deviations, which can be helpful in telling you if you there are days/times/triads
    that are not behaving as wanted.
'''))
    command.add_argument('filepath', help = ('Directory of .npz files generated from heracasa'))
    command.add_argument('working_directory', help = ('where aligned .npz files will be put.'))
    command.add_argument('--output_file',required=False, default='output_bin',metavar='O', type=str,
                         help='Output file name.')
    command.add_argument('--sigma_clip', action='store_true', required=False,
                         help='Flags triads based on a sigma clipping flagging strategy.')
    command.add_argument('--lst_start', required=True, metavar='S', type = float,
                         help='Start of LSTs, in fractional hours (3.2, 23.1 etc).')
    command.add_argument('--lst_end', required=True, metavar='R', type = float,
                         help='End of LSTs, in fractional hours (3.4, 23.5 etc).')
    command.add_argument('--date_start', required=True, metavar='R', type = int,
                         help='Start date for alignment.')
    command.add_argument('--date_end', required=True, metavar='R', type=int,
                         help='End date for alignment.')
    command.add_argument('--exclude_days', required=False, metavar='E', nargs='*', default='', help='Explicitly excluse any days in range. Specify as "2458102 2458114" etc')
    command.add_argument('--channel_number',required=False,default=1024, metavar='C', type = int,
                         help='Number of channels')
    args = command.parse_args()

    if (os.path.exists(args.working_directory) == False):
        os.makedirs(args.working_directory)
        
    # Find all .npz files. Assume they are all from heracasa...
    files = []
    for file in os.listdir(args.filepath):

        if file.endswith(".npz"):
            files.append(file)

    if args.exclude_days == '':
        parser = Julian_Parse(files,args.date_start,args.date_end,args.exclude_days)
    else:
        print(args.exclude_days)
        exclude_days_list=args.exclude_days[0].split(' ')
        parser = Julian_Parse(files,args.date_start,args.date_end,exclude_days_list)
    parser.break_up_datestamps()
    files, dates = parser.return_datestamps()

    print("Number of days: %d"%len(dates))
    # Instantiate LST_alignment class and align timestamps.
    print("Aligning LST's (use first day as reference)...")
    aligner = LST_Alignment(args.filepath,dates,files,args.lst_start,args.lst_end)
    aligner.align_timestamps()
    print("Saving Files...",end="")
    aligner.save_binned_lsts(args.working_directory+args.output_file+".npz")
    print("done")

if __name__=="__main__":
    main()
