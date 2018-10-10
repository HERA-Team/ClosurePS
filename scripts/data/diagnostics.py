import os, numpy, argparse
import casac, casa, pyuvdata
from heracasa.data.uvconv import *
from heracasa.closure import clquants

import plotly.offline as offpy
import plotly.plotly as py
from plotly.graph_objs import *
import networkx




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

    return numpy.median(closures[:,0,:,timestamp],axis=0)


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
            y = concentrations[:,scatter]
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
                      title = "Frequency (Hz)"),
                  yaxis = dict(
                      title = "|R|"))
    
    fig = dict(data=conc_plots,
               layout=layout)

    offpy.plot(fig, filename=filepath)


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
        y = concentrations
        )

    layout = dict(title = "Average Closure Concentrations over all MS Timestamps",
                  xaxis = dict(
                      title = "Frequency (Hz)"),
                  yaxis = dict(
                      title = "|R|"))
    fig = dict(data=[scatter], layout=layout)

    #py.image.save_as(fig,filename=filepath,format='png')
    offpy.plot(fig,filename=filepath,image='png')


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

    offpy.plot([scatter], filename=filepath)
                                
        

    

def process_ms_file(msname,working_dir):

    ms = casac.casac.ms()
    ms.open(msname)
    
    triangles = clquants.equiTr(msname, 14, 1)
    closures = clquants.closurePhTriads(msname, triangles)


    #med_x, med_y = plot_unit_circle(median)
    #Calculate MSE:
    mse_val_t = numpy.empty(shape=(numpy.shape(closures[:,0,0,0])))

    for j in range(numpy.shape(closures[0,0,0,:])[0]):
        median = find_median_closure(closures,j)
        other_closures = closures[:,0,:,j]
        mse_val = numpy.empty(shape=(numpy.shape(other_closures)[0]))
        for i,closure in enumerate(other_closures):
            #closure_x, closure_y = plot_unit_circle(closure)
            deviations = numpy.arctan2(numpy.sin(closure-median),numpy.cos(closure-median))
     #       print deviations
            mse_val[i] = numpy.mean(numpy.abs(deviations))
        mse_val_t = mse_val_t + mse_val
    mse_val_t = mse_val_t / numpy.shape(closures[0,0,0,:])[0]


    
    #Mask the closures to rule out any over a certain number of radians.
    mse_threshold = 3.0
    #print numpy.argwhere(mse_val_t > mse_threshold)
    resultant_closures = numpy.ma.masked_array(closures)
    resultant_closures[mse_val_t > mse_threshold] = numpy.ma.masked
    #print(numpy.shape(resultant_closures))

    ##Plot network graph:
    
    #Take first indices as antenna position:
    trigraph = triangles[:,0]
    xy = numpy.empty(shape=(len(trigraph),2))
    for i,tri in enumerate(trigraph):
        xy[i] = numpy.array(locations[tri])
        if triangles[i][0]/10 == triangles[i][1]/10:
            xy[i] += 2
    #print xy

    #Build Graph
    G = networkx.Graph()
    for i,triangle in enumerate(triangles):
        intersection, indexes = connecting_triads(triangle, triangles)
        mse_edges = calculate_mae_edges(closures, triangle, intersection, indexes, i)
        
        G.add_node(i,pos=xy[i],mse=mse_val_t[i])
        for j,edge in enumerate(indexes):
            a_coord = triangles[edge][0]
            a_xy = locations[a_coord]
            b_coord = triangles[i][0]
            b_xy = locations[b_coord]
            G.add_edge(i,edge,pos=(a_xy,b_xy),potential=mse_edges[j])

    edge_trace = Scatter(
        x=[],
        y=[],
        text=[],
        line=Line(width=0.5, color='#888'),
        hoverinfo='text',
        mode='lines')

    for edge in G.edges():
        x0, y0 = G.node[edge[0]]['pos']
        x1, y1 = G.node[edge[1]]['pos']
        #print x0, y0, x1, y1
        edge_trace['x'] += [x0, x1, None]
        edge_trace['y'] += [y0, y1, None]
        

    node_trace = Scatter(
    x=[],
    y=[],
    text=[],
    mode='markers',
    hoverinfo='text',
    marker=Marker(
        showscale=True,
        # colorscale options
        # 'Greys' | 'Greens' | 'Bluered' | 'Hot' | 'Picnic' | 'Portland' |
        # Jet' | 'RdBu' | 'Blackbody' | 'Earth' | 'Electric' | 'YIOrRd' | 'YIGnBu'
        colorscale='YIGnBu',
        reversescale=True,
        color=[],
        size=10,
        colorbar=dict(
            thickness=15,
            title='Mean Absolute Circular Error',
            xanchor='left',
            titleside='right'
        ),
        line=dict(width=2)))

    for node in G.nodes():
        x = G.node[node]['pos'][0]
        y = G.node[node]['pos'][1]
        node_trace['x'].append(x)
        node_trace['y'].append(y)
        node_trace['text'].append(G.node[node]['mse'])
        node_trace['marker']['color'].append(G.node[node]['mse'])

    fig = Figure(data=Data([edge_trace,node_trace]),
                 layout = Layout(
                     title = 'Circular MAE Values for Closures',
                     titlefont = dict(size=16),
                     showlegend = False,
                     hovermode = 'closest',
                     margin = dict(b=20,l=15,r=5,t=40)))

    offpy.plot(fig, filename=working_dir+'/graph.html')
        
    plot_concentration_measures(resultant_closures,working_dir+"/concentrations.png")
    #plot_concentration_measures_flat(resultant_closures)
    #print(numpy.shape(resultant_closures))
    plot_concentration_measures_all(resultant_closures,working_dir+"/concentrations_all.html")

#fin - Miriad File
#fout - Where UVFits and 
def mkuvfits(fin, fout):
    UV = pyuvdata.UVData()

    #Read file and rewrap phases.
    UV.read_miriad(fin, 'miriad')
    UV.phase_to_time(UV.time_array[0])
    UV.write_uvfits(fout, 'uvfits')

    #Output to .ms directory
    casa.importuvfits(fout, fout+'.ms')
    
    

def main():

    args = None
    command = argparse.ArgumentParser()
    command.add_argument('miriadfile', help = ('Miriad File Path'))
    command.add_argument('tri_size', help = ('Equilateral Triangle Size'))
    command.add_argument('working_directory', help = ('Working Directory. rds/hpc-work hopefully.'))
    args = command.parse_args()

    if (os.path.exists(args.working_directory) == False):
        os.makedirs(args.working_directory)

    path = args.miriadfile.strip("/")
    basename = os.path.basename(path)
    
    #Convert miriad file to ms file
    mkuvfits(args.miriadfile, args.working_directory+basename)
    process_ms_file(args.working_directory+basename+".ms",args.working_directory)
    



    
if __name__=="__main__":
    main()




