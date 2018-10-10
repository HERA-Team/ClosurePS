import os, glob, numpy, argparse
from astropy.time import Time

import multiprocessing


import plotly.offline as offpy
import plotly.plotly as py
from plotly.graph_objs import *
import networkx

# Number of cores on node. CSD3 has 32 per node.
# This is IO limited so filling up all the cores will reduce performance.
NUMBER_OF_CORES = 8

#Triads that we are processing. 
triangles = numpy.array([[  0,   1,  12],
             [  1,   2,  13],
             [ 11,  12,  24],
             [ 12,  13,  25],
             [ 13,  14,  26],
             [ 23,  24,  37],
             [ 24,  25,  38],
             [ 25,  26,  39],
             [ 26,  27,  40],
             [ 36,  37,  51],
             [ 37,  38,  52],
             [ 38,  39,  53],
             [ 40,  41,  55],
             [ 51,  52,  67],
             [ 52,  53,  68],
             [ 53,  54,  69],
             [ 54,  55,  70],
             [ 65,  66,  82],
             [ 66,  67,  83],
             [ 67,  68,  84],
             [ 68,  69,  85],
             [ 70,  71,  87],
             [120, 121, 140],
             [121, 122, 141],
             [122, 123, 142],
             [123, 124, 143]])

#Fake co-ordinates for every antenna in HERA-47, for plotting.
locations = dict([(0,(-10,-50)),
                  (1,(0,-50)),
                  (2,(10,-50)),
                  (11,(-15,-40)),
                  (12,(-5,-40)),
                  (13,(5,-40)),
                  (14,(15,-40)),
                  (23,(-20,-30)),
                  (24,(-10,-30)),
                  (25,(0,-30)),
                  (26,(10,-30)),
                  (27,(20,-30)),
                  (36,(-25,-20)),
                  (37,(-15,-20)),
                  (38,(-5,-20)),
                  (39,(5,-20)),
                  (40,(15,-20)),
                  (41,(25,-20)),
                  (50,(-30,-10)),
                  (51,(-20,-10)),
                  (52,(-10,-10)),
                  (53,(0,-10)),
                  (54,(10,-10)),
                  (55,(20,-10)),
                  (65,(-35,0)),
                  (66,(-25,0)),
                  (67,(-15,0)),
                  (68,(-5,0)),
                  (69,(5,0)),
                  (70,(15,0)),
                  (71,(25,0)),
                  (82,(-30,10)),
                  (83,(-20,10)),
                  (84,(-10,10)),
                  (85,(0,10)),
                  (86,(10,10)),
                  (87,(20,10)),
                  (88,(30,10)),
                  (98,(-45,20)),
                  (120,(-10,30)),
                  (121,(0,30)),
                  (122,(10,30)),
                  (123,(20,30)),
                  (124,(30,30)),
                  (141,(5,40)),
                  (142,(15,40)),
                  (143,(25,40))
])

msname = "/home/jck42/casascripts/hera_alt.ms/"

def plot_unit_circle(closures_chan):

    xc = numpy.cos(closures_chan)
    yc = numpy.sin(closures_chan)

    return xc, yc

def find_median_closure(closures,timestamp):

    return numpy.ma.median(closures[:,0,:,timestamp],axis=0)


#Find intersections between triads and other sets.
#This finds the edges between the nodes, allowing us to construct our graph.
def connecting_triads(triangle, triangles):

    intersection = []
    indexes = []
    for  i,tri  in enumerate(triangles):
        inter = numpy.intersect1d(triangle,tri)
        if inter.size > 0 and inter.size<3:
            intersection.append(tri)
            indexes.append(i)

    return numpy.array(intersection), numpy.array(indexes)

def calculate_mae_edges(closures, triangle, intersections, index_closures, index_centre):    

    closure_triangle = closures[index_centre,0,:,0]
    mse = []
    for i in index_closures:
        closures_intersect = closures[i,0,:,0]
        mse.append(numpy.mean(numpy.abs(numpy.arctan2(numpy.sin(closures_intersect-closure_triangle),
                                                      numpy.cos(closures_intersect-closure_triangle)))))
#        mse.append(calculate_mse(closure_triangle, closures_intersect, 1024))


    return numpy.array(mse)

def find_uv_of_triangles(triangles):

    pass

def conc_measure(angular_data):
    
    x, y = plot_unit_circle(angular_data)
    concentration = numpy.sqrt(numpy.square(numpy.mean(x)) + numpy.square(numpy.mean(y)))
    return concentration


def circular_average(closures):

    x_var = numpy.cos(closures)
    y_var = numpy.sin(closures)

    average_angle = numpy.arctan2(numpy.mean(y_var),numpy.mean(x_var))
    return average_angle

#Plots all concentrations with a slider bar to navigate.
def plot_concentration_measures_all(closures,filepath):

    concentrations = numpy.empty(shape=(numpy.shape(closures[0,0,:,:])))
    cbw = 1.0e8 / 1024
    for i in range(numpy.shape(closures)[2]):
        for j in range(numpy.shape(closures)[3]):
            concentrations[i,j] = conc_measure(closures[:,0,i,j])
    conc_plots = []

    print(concentrations[:,3])

    for scatter in range(numpy.shape(closures)[3]):
        
        scatter = Scatter(
            visible = False,
            x = 1.0e8 + cbw * numpy.arange(0, numpy.shape(closures)[2]),
            y = concentrations[:,scatter],
            mode = 'markers'
        )
        conc_plots.append(scatter)

    timestamps = []
    # print concentrations[:,0]
    for i in range(numpy.shape(closures)[3]):
        timestamp = dict(
            method = 'restyle',
            args = ['visible', [False] * numpy.shape(concentrations)[1]],
            label = i
        )
        timestamp['args'][1][i] = True # Toggle i'th trace to "visible"
        timestamps.append(timestamp)

    sliders = [dict(
        active = 100,
        currentvalue = {"prefix": "Closure Timestamp: "},
        #pad = {"t": 50},
        steps = timestamps
    )]

    layout = dict(sliders=sliders,
                  title = "Closure Concentrations per MS Timestamp",
                  xaxis = dict(
                      title = "Frequency (Hz)",
                      range = [1.0e8, 2.0e8]),
                  yaxis = dict(
                      title = "|R|",
                      range = [0 ,1]),
    )
    
    fig = dict(data=conc_plots,
               layout=layout)

    offpy.plot(fig, filename=filepath, auto_open=False)


#Plots concentration vector for all triads vs channel number. Concentration vector
#is calculated for each timestamp. These individual vectors are averaged together at the end.
def plot_concentration_measures(closures, filepath):

    #py.sign_in('JamesKent','qGRzuezlobl2BLMlp4iZ')
    #Average over all timestamps for each channel
    concentrations = numpy.empty(shape=(numpy.shape(closures)[2]))
    cbw = 1.0e8 / 1024
    for i in range(numpy.shape(closures)[2]):
        conct = 0.0
        for j in range(numpy.shape(closures)[3]):
            conct += conc_measure(closures[:,0,i,j])
        #print conct
        conct = conct / numpy.shape(closures)[3]
        #print conct
        concentrations[i] = conct

    scatter = Scatter(
        x = 1.0e8 + numpy.arange(0,numpy.shape(closures)[2]) * cbw,
        y = concentrations,
        mode = 'markers'
        )

    layout = dict(title = "Average Closure Concentrations over all MS Timestamps",
                  xaxis = dict(
                      title = "Frequency (Hz)",
                      range = [1.0e8, 2.0e8]),
                  yaxis = dict(
                      title = "|R|",
                      range = [0, 1]))
    fig = dict(data=[scatter], layout=layout)

    #py.image.save_as(fig,filename=filepath,format='png')
#    offpy.plot(fig,filename=filepath,image='png')
    offpy.plot(fig,filename=filepath, auto_open=False)


#Plots concentration vector across all timestamps of an ms file (likely not that useful)    
def plot_concentration_measures_flat(closures, filepath):

    concentrations = numpy.empty(shape=numpy.shape(closures)[2])
    for i in range(numpy.shape(closures)[2]):

        x_tot = numpy.empty(shape=(numpy.shape(closures)[0]))
        y_tot = numpy.empty(shape=(numpy.shape(closures)[0]))
        for j in range(numpy.shape(closures)[3]):
            x,y = plot_unit_circle(closures[:,0,i,j])
            x_tot += x
            y_tot += y

        x_m = numpy.sum(x_tot)/(numpy.shape(closures)[2] * numpy.shape(closures)[3])
        y_m = numpy.sum(y_tot)/(numpy.shape(closures)[2] * numpy.shape(closures)[3])
        concentrations[i] = numpy.sqrt(numpy.square(x_m)+numpy.square(y_m))

    scatter = Scatter(
        x = numpy.arange(0, numpy.shape(closures)[2]),
        y = concentrations)

    offpy.plot([scatter], filename=filepath, auto_open=False)
                                
        

def plot_closures_channel(closures, filepath, timestamp = 0):

    cbw = 1.0e8 / 1024
    plots = []
    for scatter in range(numpy.shape(closures)[0]):

        scatter = Scatter(
            visible = True,
            x = 1.0e8 + cbw * numpy.arange(0, numpy.shape(closures)[2]),
            y = closures[scatter,0,:,timestamp],
            mode = 'markers'
            )
        plots.append(scatter)

    layout = dict(title = "All Closures at timestamp "+str(timestamp),
                  xaxis = dict (
                      title = "Frequency (Hz)",
                      range = [1.0e8, 2.0e8]),
                  yaxis = dict (
                      title = "Phase",
                      range = [-numpy.pi,numpy.pi])
    )
    fig = dict(data=plots,
               layout=layout)

    offpy.plot(fig, filename=filepath, image='png',auto_open=False)


def plot_closures_average_all(closures, filepath):

    
    cbw = 1.0e8 / 1024
    closure_average = numpy.empty(shape=numpy.shape(closures[0,0,:,:]))
    for i in range(numpy.shape(closures)[2]):
        for j in range(numpy.shape(closures)[3]):
            closure_average[i,j] = circular_average(closures[:,0,i,j])

    concentrations = numpy.empty(shape=numpy.shape(closures[0,0,:,:]))
    for i in range(numpy.shape(closures)[2]):
        for j in range(numpy.shape(closures)[3]):
            concentrations[i,j] = conc_measure(closures[:,0,i,j])
    std_devs = numpy.sqrt(-2 * numpy.log(concentrations))
    
    
    plots = []

    for scatter_i in range(numpy.shape(closures)[3]):
        scatter = Scatter(
            visible = False,
            name = 'Circular Average',
            x = 1.0e8 + cbw * numpy.arange(0, numpy.shape(closures)[2]),
            y = closure_average[:,scatter_i],
            mode = 'markers'
        )
        plots.append(scatter)

            
        variances = Scatter(
            visible = False,
            name = 'Circular Std. Dev.',
            x = 1.0e8 + cbw * numpy.arange(0, numpy.shape(closures)[2]),
            y = std_devs[:,scatter_i],
            mode = 'markers'
        )
        plots.append(variances)


    timestamps = []
    for i in range(numpy.shape(closures)[3]):
        timestamp = dict(
            method = 'restyle',
            args = ['visible', [False] * 2 * numpy.shape(concentrations)[1]],
            label = i
        )
        timestamp['args'][1][2 * i] = True # Toggle i'th trace to "visible"
        timestamp['args'][1][2 * i + 1] = True
        timestamps.append(timestamp)

    sliders = [dict(
        active = 100,
        currentvalue = {"prefix": "Closure Timestamp: "},
        #pad = {"t": 50},
        steps = timestamps
    )]
    
        
    layout = dict(sliders=sliders,
                  title = "Average Closures at Timestamp",
                  xaxis = dict(
                      title = "Frequency(Hz)",
                      range=[1.0e8, 2.0e8]),
                  yaxis = dict(
                      title = "Phase(Radians)",
                      range =[-numpy.pi,numpy.pi]))

    fig = dict(data = plots,
               layout = layout)

    offpy.plot(fig, filename=filepath, auto_open=False)


def plot_closures_average(closures, filepath, timestamp = 0):

    cbw = 1.0e8 / 1024
    closure_ts = closures[:,0,:,timestamp]
    
    closure_average = numpy.empty(shape=numpy.shape(closures)[2])
    for i in range(numpy.shape(closures)[2]):
        closure_average[i] = circular_average(closures[:,0,i,timestamp])

    concentrations = numpy.empty(shape=numpy.shape(closures)[2])
    for i in range(numpy.shape(closures)[2]):
        concentrations[i] = conc_measure(closures[:,0,i,timestamp])

    std_devs = numpy.sqrt(-2 * numpy.log(concentrations))
    
    print(numpy.shape(closure_average))
        

    plots = []
    scatter = Scatter(
        visible = True,
        name = 'Circular Average',
        x = 1.0e8 + cbw * numpy.arange(0, numpy.shape(closures)[2]),
        y = closure_average,
        mode = 'markers'
        )
    plots.append(scatter)

    variances = Scatter(
        visible = True,
        name = 'Circular Std. Dev.',
        x = 1.0e8 + cbw * numpy.arange(0, numpy.shape(closures)[2]),
        y = std_devs,
        mode = 'markers'
        )
    plots.append(variances)

    layout = dict(title = "Average Closures at Timestamp",
                  xaxis = dict(
                      title = "Frequency(Hz)"),
                  yaxis = dict(
                      title = "Phase(Radians)"))

    fig = dict(data = plots,
               layout = layout)

    offpy.plot(fig, filename=filepath, image='png',auto_open=False)
    

def process_npz_file(filename,working_dir):

    # Load closures from npz file.
    with numpy.load(filename) as data:
        closures = data['phase']
        flags = data['flags']

    
    ## MASKING:
    #Mask out any channels that have been flagged.
    #Also mask out channels above 180MHz and below 115 MHz.
    #Closures go a bit bananas outside this range.
    closures = numpy.ma.masked_array(closures)
    closures[flags == True] = numpy.ma.masked
    closures[:,:,:150,:] = numpy.ma.masked
    closures[:,:,819:,:] = numpy.ma.masked
    #Unique UID for each output file.
    file_uid = os.path.basename(filename)

    
    
    #Calculate MSE:
    mse_val_t = numpy.empty(shape=(numpy.shape(closures[:,0,0,0])))

    for j in range(numpy.shape(closures[0,0,0,:])[0]):
        median = find_median_closure(closures,j)
        other_closures = closures[:,0,:,j]
        mse_val = numpy.empty(shape=(numpy.shape(other_closures)[0]))
        for i,closure in enumerate(other_closures):
            deviations = numpy.arctan2(numpy.sin(closure-median),numpy.cos(closure-median))
            mse_val[i] = numpy.ma.mean(numpy.abs(deviations))
            print(mse_val[i])
            
        mse_val_t = mse_val_t + mse_val
    mse_val_t = mse_val_t / numpy.shape(closures[0,0,0,:])[0]
    

    
    #Mask the closures to rule out any over a certain number of radians.
    mse_threshold = 3.0
    resultant_closures = numpy.ma.masked_array(closures)
    resultant_closures[mse_val_t > mse_threshold] = numpy.ma.masked
    ##Plot network graph:
    
    # #Take first indices as antenna position:
    # trigraph = triangles[:,0]
    # xy = numpy.empty(shape=(len(trigraph),2))
    # for i,tri in enumerate(trigraph):
    #     xy[i] = numpy.array(locations[tri])
    #     if triangles[i][0]/10 == triangles[i][1]/10:
    #         xy[i] += 2
    # #print xy

    # #Build Graph
    # G = networkx.Graph()
    # for i,triangle in enumerate(triangles):
    #     intersection, indexes = connecting_triads(triangle, triangles)
    #     #print(mse_val_t[i])
    #     G.add_node(i,pos=xy[i],mse=mse_val_t[i])
    #     for j,edge in enumerate(indexes):
    #         a_coord = triangles[edge][0]
    #         a_xy = locations[a_coord]
    #         b_coord = triangles[i][0]
    #         b_xy = locations[b_coord]
    #         G.add_edge(i,edge,pos=(a_xy,b_xy))

    # edge_trace = Scatter(
    #     x=[],
    #     y=[],
    #     text=[],
    #     line=Line(width=0.5, color='#888'),
    #     hoverinfo='text',
    #     mode='lines')

    # for edge in G.edges():
    #     x0, y0 = G.node[edge[0]]['pos']
    #     x1, y1 = G.node[edge[1]]['pos']
    #     #print x0, y0, x1, y1
    #     edge_trace['x'] += [x0, x1, None]
    #     edge_trace['y'] += [y0, y1, None]
        

    # node_trace = Scatter(
    # x=[],
    # y=[],
    # text=[],
    # mode='markers',
    # hoverinfo='text',
    # marker=Marker(
    #     showscale=True,
    #     # colorscale options
    #     # 'Greys' | 'Greens' | 'Bluered' | 'Hot' | 'Picnic' | 'Portland' |
    #     # Jet' | 'RdBu' | 'Blackbody' | 'Earth' | 'Electric' | 'YIOrRd' | 'YIGnBu'
    #     colorscale='YIGnBu',
    #     reversescale=True,
    #     color=[],
    #     size=10,
    #     colorbar=dict(
    #         thickness=15,
    #         title='Mean Absolute Circular Error',
    #         xanchor='left',
    #         titleside='right'
    #     ),
    #     line=dict(width=2)))

    # for node in G.nodes():
    #     x = G.node[node]['pos'][0]
    #     y = G.node[node]['pos'][1]
    #     node_trace['x'].append(x)
    #     node_trace['y'].append(y)
    #     node_trace['text'].append(G.node[node]['mse'])
    #     node_trace['marker']['color'].append(G.node[node]['mse'])

    # fig = Figure(data=Data([edge_trace,node_trace]),
    #              layout = Layout(
    #                  title = 'Circular MAE Values for Closures',
    #                  titlefont = dict(size=16),
    #                  showlegend = False,
    #                  hovermode = 'closest',
    #                  margin = dict(b=20,l=15,r=5,t=40)))

#    offpy.plot(fig, filename=working_dir+"/"+file_uid+'_graph.html',auto_open=False)
    plot_closures_channel(resultant_closures,working_dir+"/"+file_uid+"_closures.html")
    plot_closures_average_all(resultant_closures,working_dir+"/"+file_uid+"_closure_average.html")       
    plot_concentration_measures(resultant_closures,working_dir+"/"+file_uid+"_concentrations.html")
    plot_concentration_measures_all(resultant_closures,working_dir+"/"+file_uid+"_concentrations_all.html")
    

def multiprocess_wrapper(filepaths, main_filepath, working_directory):

    for file in filepaths:
        process_npz_file(main_filepath+file,working_directory)


    
def main():

    args = None
    command = argparse.ArgumentParser()
    command.add_argument('filepath', help = ('.npz directory'))
    command.add_argument('working_directory', help = ('Working Directory. Where files will be put.'))
    args = command.parse_args()

    if (os.path.exists(args.working_directory) == False):
        os.makedirs(args.working_directory)

    files = os.listdir(args.filepath)
    npz_files =[]
    for file in files:
        if file.endswith(".npz"):
            npz_files.append(file)


    #Process files in seperate processes.
            
    files_per_core = len(npz_files) / NUMBER_OF_CORES
    split_files = [npz_files[i:i+files_per_core] for i in range(0, len(npz_files), files_per_core)]
    remainder_files = npz_files[files_per_core * NUMBER_OF_CORES:len(npz_files)]

    jobs = []

    for i, list_slice in enumerate(split_files):
        print(list_slice)
        j = multiprocessing.Process(target=multiprocess_wrapper, args=(list_slice,args.filepath,args.working_directory))
        jobs.append(j)

    for job in jobs:
        job.start()
                            

    # Process remainder in host process.
    for file in remainder_files:
        process_npz_file(args.filepath+file,args.working_directory)
    
    
    #for file in os.listdir(args.filepath):
    #    if file.endswith(".npz"):
    #        process_npz_file(args.filepath+file,args.working_directory)


    
if __name__=="__main__":
    main()




