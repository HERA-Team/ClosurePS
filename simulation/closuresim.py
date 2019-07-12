import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D 

# For co-ordinate transforms

def radec_to_hadec(ra,dec,lst):    
    return lst-ra,dec

def hadec_to_radec(ha,dec,lst):
    return lst-ha,dec

def hadec_to_altaz(ha,dec,lat):
    assert(numpy.abs(lat)<numpy.pi/2)
    assert(numpy.any(dec > -90) and numpy.any(dec < 90))

    alt = numpy.arcsin(numpy.sin(dec) * numpy.sin(lat) + numpy.cos(dec) * numpy.cos(lat) * numpy.cos(ha))
    a = numpy.arccos((numpy.sin(dec) - numpy.sin(lat) * numpy.sin(alt))/(numpy.cos(lat) * numpy.cos(alt)))
    
    tr = numpy.sin(ha) > 0
    a[tr] = (2 * numpy.pi) - a[tr] 
    
    
    return alt,a

def altaz_to_hadec(alt,az,lat):
    dec = numpy.arcsin(numpy.sin(alt) * numpy.sin(lat) + numpy.cos(alt) * numpy.cos(lat) * numpy.cos(az))
    ha = numpy.arccos((numpy.sin(alt) - numpy.sin(lat) * numpy.sin(dec))/(numpy.cos(lat) * numpy.cos(dec)))    
    return ha,dec

def altaz_to_dircos(alt,az):
    phi = numpy.pi/2 - az
    l = numpy.cos(alt) * numpy.cos(phi)
    m = numpy.cos(alt) * numpy.sin(phi)
    n = numpy.sqrt(1.0 - l*l - m*m) - 1.0
    return l,m,n
def dircos_to_altaz(l,m):
    
    n = numpy.sqrt(1.0 - l*l - m*m) - 1.0
    alt = numpy.pi/2 - numpy.arccos(n)
    az = numpy.pi/2 - numpy.arctan2(m,l)
    return alt,az


# For simulating Electric Fields
def simulate_electric_fields(locations_uvw,points):
    elec_ss = numpy.zeros(shape=(locations_uvw.shape[0]),dtype=numpy.complex128)
    
    for i in numpy.arange(locations_uvw.shape[0]):
        for point in points:
            l = point[0]
            m = point[1]
            n = numpy.sqrt(1.0-l**2-m**2) - 1.0
            amp = point[2]
            u = locations_uvw[i,0]
            v = locations_uvw[i,1]
            w = locations_uvw[i,2]
            elec_ss[i] += amp*numpy.exp(-2j*numpy.pi*(u*l+v*m+w*n))
    return elec_ss

# 7 Antennas
def generate_hera_layout_simple():
    layout_vec = numpy.zeros(shape=(14,3))
    
    
    layout_vec[0] = numpy.asarray([-7,24.2486,0])
    layout_vec[1] = numpy.asarray([7,24.2486,0])
    layout_vec[2] = numpy.asarray([-14,12.1243,0])
    layout_vec[3] = numpy.asarray([0,12.1243,0])
    layout_vec[4] = numpy.asarray([14,12.1243,0])
    layout_vec[5] = numpy.asarray([-21,0,0])
    layout_vec[6] = numpy.asarray([-7,0,0])
    layout_vec[7] = numpy.asarray([7,0,0])
    layout_vec[8] = numpy.asarray([21,0,0])
    layout_vec[9] = numpy.asarray([-14,-12.1243,0])
    layout_vec[10] = numpy.asarray([0,-12.1243,0])
    layout_vec[11] = numpy.asarray([14,-12.1243,0])
    layout_vec[12] = numpy.asarray([-7,-24.2486,0])
    layout_vec[13] = numpy.asarray([7,-24.2486,0])

    
    
    return layout_vec

from scipy.stats import norm

def generate_gaussian_beam():

    x = numpy.linspace(-1,1,2000)
    return x,norm.pdf(x, loc=0, scale=0.2)+ 1j*norm.pdf(x, loc=0,scale=0.2)#*norm.cdf(2*x,loc=0,scale=0.2) + 1j*norm.pdf(x, loc=0,scale=0.2)#*norm.cdf(2*x,loc=0,scale=0.2)

def find_beam_at_point(beam,point):
    beam_shape = beam.shape[0]
    
    pl = point[0]
    pm = point[1]
    
    step = 2/beam_shape
    
    bl_ind = int((pl + 1)/step)
    bm_ind = int((pm + 1)/step)
    
    # Multiple real/imaginary seperately. 
    
    re_be = beam[bl_ind].real * beam[bm_ind].real
    im_be = beam[bl_ind].imag * beam[bm_ind].imag    
    
    return re_be + 1j * im_be

    


# For simulating Electric Fields
def simulate_electric_fields_gaussianbeam(locations_uvw,points,beam):
    elec_ss = numpy.zeros(shape=(locations_uvw.shape[0]),dtype=numpy.complex128)
    phase_source = numpy.zeros(shape=(points.shape[0]),dtype=numpy.complex128)
    print(phase_source.shape)
    for i in numpy.arange(points.shape[0]):
        phase_source[i] = numpy.exp(2*numpy.pi * numpy.random.uniform(0,1))
    
    for i in numpy.arange(locations_uvw.shape[0]):
        #print(locations_uvw[i,:])
        #numpy.random.seed(12345)
        #antenna_noise = 500*(numpy.random.normal(scale=4000.0) + 1j * numpy.random.normal(scale=4000.0))
        #antenna_phase_noise = numpy.exp(1j*numpy.random.uniform(-0.1,0.1))
        #antenna_noise = 0
        
        for j,point in enumerate(points):
            #print(j)
            l = point[0]
            m = point[1]
            n = numpy.sqrt(1.0-l**2-m**2) - 1.0
            u = locations_uvw[i,0]
            v = locations_uvw[i,1]
            w = locations_uvw[i,2]
            elec_ss[i] += point[2]*numpy.exp(-2j*numpy.pi*(u*l+v*m+w*n))  * find_beam_at_point(beam,point)# * numpy.exp(2*numpy.pi * numpy.random.uniform(0,1)) #* numpy.random.uniform(0.9,1.1)
        #elec_ss[i] += antenna_noise #* antenna_phase_noise
    return elec_ss

def epic_image_dft(theta,lam,efield,locations):
    grid_size = int(theta*lam)
    image = numpy.zeros(shape=(grid_size,grid_size),dtype=numpy.complex128)
    
    for lg in numpy.arange(grid_size):
        
        lv = (1.0 / (0.5*grid_size)) * (lg - grid_size/2)
        print(lv)
        for mg in numpy.arange(grid_size):
            mv = (1.0/ (0.5*grid_size)) * (mg - grid_size/2)
            if (lv**2 + mv**2 > 1.0):
                continue
                
            else:
                nv = 1.0 - numpy.sqrt(1.0 - lv**2 - mv**2)
            #print(nv)
            for j,efieldd in enumerate(efield):
                u = locations[j,0]
                v = locations[j,1]
                w = locations[j,2]
                #w = 0
                image[lg,mg] += efieldd * numpy.exp(2j*numpy.pi*(u*lv + v*mv + w*nv))
    return image / (grid_size**2)

# For simulating Electric Fields
def simulate_electric_fields_gaussianbeam(locations_uvw,points,beam):
    elec_ss = numpy.zeros(shape=(locations_uvw.shape[0]),dtype=numpy.complex128)
    phase_source = numpy.zeros(shape=(points.shape[0]),dtype=numpy.complex128)
    print(phase_source.shape)
    for i in numpy.arange(points.shape[0]):
        phase_source[i] = numpy.exp(2*numpy.pi * numpy.random.uniform(0,1))
    
    for i in numpy.arange(locations_uvw.shape[0]):
        #print(locations_uvw[i,:])
        #numpy.random.seed(12345)
        #antenna_noise = 500*(numpy.random.normal(scale=4000.0) + 1j * numpy.random.normal(scale=4000.0))
        #antenna_phase_noise = numpy.exp(1j*numpy.random.uniform(-0.1,0.1))
        #antenna_noise = 0
        
        for j,point in enumerate(points):
            #print(j)
            l = point[0]
            m = point[1]
            n = numpy.sqrt(1.0-l**2-m**2) - 1.0
            u = locations_uvw[i,0]
            v = locations_uvw[i,1]
            w = locations_uvw[i,2]
            elec_ss[i] += point[2]*numpy.exp(-2j*numpy.pi*(u*l+v*m+w*n))  * find_beam_at_point(beam,point)# * numpy.exp(2*numpy.pi * numpy.random.uniform(0,1)) #* numpy.random.uniform(0.9,1.1)
        #elec_ss[i] += antenna_noise #* antenna_phase_noise
    return elec_ss


#### VISIBILITIES #####

def generate_baseline_vectors(locations):
    
    loc_no = locations.shape[0]
    baseline_vectors = numpy.zeros(shape=(loc_no, loc_no,3))
    
    for i in numpy.arange(loc_no):
        for j in numpy.arange(loc_no):
            baseline_vectors[i,j,0] = (locations[i,0] - locations[j,0])#/2
            baseline_vectors[i,j,1] = (locations[i,1] - locations[j,1])#/2
            baseline_vectors[i,j,2] = (locations[i,2] - locations[j,2])#/2
            
    return baseline_vectors


def plot_baselines(baselines):
    
    us = baselines[:,:,0]
    vs = baselines[:,:,1]
    ws = baselines[:,:,2]
    ax = plt.subplot(121, projection='3d')
    ax.scatter(us,vs,ws, color='red')
    ax.set_title('Electric Field Locations')
    ax.set_xlabel('u'); ax.set_ylabel('v'); ax.set_zlabel('w')
    plt.show()
    
def generate_gaussian_power_beam(mean,sigma):
    x = numpy.linspace(-1,1,2000)
    pdf = numpy.exp(-(((x-mean)**2)/(2*sigma**2)))#*(1/numpy.sqrt(2*numpy.pi*sigma**2))
    return x,pdf

def generate_gaussian_power_beams(locs,mean,lb,ub,lbp,ubp):
    
    # Locs = antenna locs in u/v/w (wavelengths)
    # Mean = Mean loc of antenna beam. Zero will suffice if centred on zenith in direction cosine co-ordinates
    # lb = Lower bound of sigmas for real part of beam
    # ub = Upper bound ""
    # lbp = Lower bound of sigmas for imaginary party of beam
    # ubp = Upper bound "" 
    
    
    assert(locs.shape[0] == mean.shape[0]) #== sigmas.shape[0])
    beams_matrix = numpy.zeros(shape=(locs.shape[0],2000),dtype=numpy.complex128)
    x = numpy.linspace(-1,1,2000)
    
    for i in numpy.arange(locs.shape[0]):
        sigma = numpy.random.uniform(lb,ub)
        sigma_imag = numpy.random.uniform(lbp,ubp)
        beams_matrix[i,:] = numpy.exp(-(((x-mean[i])**2)/(2*sigma**2)))# + 1j*numpy.exp(-(((x-mean[i])**2)/(2*sigma_imag**2)))
        
        # Normalise to unit amplitude
        beam_peak = numpy.max(numpy.abs(beams_matrix[i,:]))
        beams_matrix[i,:] /= beam_peak
        
    mean_beam = numpy.mean(beams_matrix,axis=0) # Beam error
    beam_error = numpy.abs(mean_beam - beams_matrix)
    av_beam_error = numpy.mean(beam_error)
    print(av_beam_error)
    

    
    return x,beams_matrix,av_beam_error
    
def simulate_visibilities(layout,points,beam):
    
    vis_matrix = numpy.zeros(shape=(layout.shape[0],layout.shape[0]),dtype=numpy.complex128)
    
    for i in numpy.arange(layout.shape[0]):
        for j in numpy.arange(layout.shape[0]):        
            
            ud = layout[i,0] - layout[j,0]
            vd = layout[i,1] - layout[j,1]
            wd = layout[i,2] - layout[j,2]
            
            
            
            for p in numpy.arange(points.shape[0]):
            
                lp = points[p,0]
                mp = points[p,1]
                np = numpy.sqrt(1.0 - lp**2+mp**2) - 1.0
                
                b1 = find_beam_at_point(beam,points[p,:]) 
                #print(b1)
                b2 = find_beam_at_point(beam,points[p,:])
                bap = b1 * numpy.conj(b2) # Each antenna has same beam so square
                
                phase_factor = numpy.exp(-2j * numpy.pi * (ud*lp + vd*mp + wd*np))
                point_amp = points[p,2]
                vis_matrix[i,j] += point_amp * bap * phase_factor
                
    return vis_matrix


def simulate_visibilities_beams(layout,points,beams):
    
    vis_matrix = numpy.zeros(shape=(layout.shape[0],layout.shape[0]),dtype=numpy.complex128)
    
    for i in numpy.arange(layout.shape[0]):
        for j in numpy.arange(layout.shape[0]):        
            
            ud = layout[i,0] - layout[j,0]
            vd = layout[i,1] - layout[j,1]
            wd = layout[i,2] - layout[j,2]
            
            
            
            for p in numpy.arange(points.shape[0]):
            
                lp = points[p,0]
                mp = points[p,1]
                np = numpy.sqrt(1.0 - lp**2+mp**2) - 1.0
                b1 = find_beam_at_point(beams[i,:],points[p,:])
                
                b2 = find_beam_at_point(beams[j,:],points[p,:])
                bap = b1 * numpy.conj(b2) # Each antenna has same beam so square
                
                phase_factor = numpy.exp(-2j * numpy.pi * (ud*lp + vd*mp + wd*np))
                point_amp = points[p,2]
                vis_matrix[i,j] += point_amp * bap * phase_factor
                
    return vis_matrix      



def simulate_visibilities_sky(layout,sky,beam):
    
    vis_matrix = numpy.zeros(shape=(layout.shape[0],layout.shape[0]),dtype=numpy.complex128)
    grid_size = sky.shape[0]
    lm_step = 2/grid_size
    
    
    for i in numpy.arange(layout.shape[0]):
        for j in numpy.arange(layout.shape[0]):        
            
            ud = layout[i,0] - layout[j,0]
            vd = layout[i,1] - layout[j,1]
            wd = layout[i,2] - layout[j,2]            
            
            for lc in numpy.arange(grid_size):
                lp = (lc-grid_size//2) * lm_step
                for mc in numpy.arange(grid_size):
                    mp = (mc-grid_size//2) * lm_step
                    if(lp**2 + mp**2 > 1):
                        continue
                    np = numpy.sqrt(1.0-lp**2-mp**2) - 1.0
                    b1 = find_beam_at_point(beam,[lp,mp,0])
                    b2 = find_beam_at_point(beam,[lp,mp,0])
                    bap = b1 * numpy.conj(b2) # Each antenna has same beam so square
                
                    phase_factor = numpy.exp(-2j * numpy.pi * (ud*lp + vd*mp + wd*np))
                    point_amp = sky[mc,lc]
                    vis_matrix[i,j] += point_amp * bap * phase_factor
                
    return vis_matrix/(grid_size**2)


def simulate_visibilities_beams_sky(layout,sky,beams):
    
    vis_matrix = numpy.zeros(shape=(layout.shape[0],layout.shape[0]),dtype=numpy.complex128)
    grid_size = sky.shape[0]
    lm_step = 2/grid_size
    
    for i in numpy.arange(layout.shape[0]):
        for j in numpy.arange(layout.shape[0]):        
            
            ud = layout[i,0] - layout[j,0]
            vd = layout[i,1] - layout[j,1]
            wd = layout[i,2] - layout[j,2]
            lm_step = 2/grid_size       
            
            for lc in numpy.arange(grid_size):
                lp = (lc-grid_size//2) * lm_step
                for mc in numpy.arange(grid_size):
                    mp = (mc-grid_size//2) * lm_step
                    if(lp**2 + mp**2 > 1):
                        continue
                    np = numpy.sqrt(1.0-lp**2-mp**2) - 1.0
                    b1 = find_beam_at_point(beams[i,:],[lp,mp,0])
                    b2 = find_beam_at_point(beams[j,:],[lp,mp,0])
                    bap = b1 * numpy.conj(b2) # Each antenna has same beam so square
                
                    phase_factor = numpy.exp(-2j * numpy.pi * (ud*lp + vd*mp + wd*np))
                    point_amp = sky[mc,lc]
                    vis_matrix[i,j] += point_amp * bap * phase_factor
                
    return vis_matrix/(grid_size**2)   

    
#### Closure Phases #### 

def generate_closure_triads(visb,antenna_list):
    
    ant_list = numpy.asarray(antenna_list)
    triads = numpy.zeros(shape=ant_list.shape[0])
    for i,triadconfig in enumerate(ant_list):
        vis1l = triadconfig[0]
        vis2l = triadconfig[1]
        vis3l = triadconfig[2]
        #print(triadconfig)
        vis1 = numpy.angle(visb[vis1l[0],vis1l[1]])
        vis2 = numpy.angle(visb[vis2l[0],vis2l[1]])
        vis3 = numpy.angle(visb[vis3l[0],vis3l[1]])
        print(vis1,vis2,vis3)
        triads[i] = vis1 + vis2 + vis3
        
    return triads

def generate_closure_triads_tp(visb,antenna_list):
    
    ant_list = numpy.asarray(antenna_list)
    triads = numpy.zeros(shape=ant_list.shape[0])
    for i,triadconfig in enumerate(ant_list):
        vis1l = triadconfig[0]
        vis2l = triadconfig[1]
        vis3l = triadconfig[2]
        #print(triadconfig)
        triple_product = numpy.angle(visb[vis1l[0],vis1l[1]] * visb[vis2l[0],vis2l[1]] * visb[vis3l[0],vis3l[1]])
        vis1 = numpy.arctan2(numpy.sin(triple_product),numpy.cos(triple_product))
        #vis2 = numpy.angle(visb[vis2l[0],vis2l[1]])
        #vis3 = numpy.angle(visb[vis3l[0],vis3l[1]])
        #print(vis1,vis2,vis3)
        triads[i] = vis1# + vis2 + vis3
        
    return triads

def calculate_closure_standard_deviation(closures):
    xc = numpy.cos(closures)
    yc = numpy.sin(closures)
    
    conc = numpy.sqrt(numpy.square(numpy.mean(xc)) + numpy.square(numpy.mean(yc)))
    return numpy.sqrt(-2 * numpy.log(conc))

def generate_beam_per_antenna(locs,lb,ub):
    beams_matrix = numpy.zeros(shape=(locs.shape[0],2000),dtype=numpy.complex128)
    
    x = numpy.linspace(-1,1,2000)
    
    for i in numpy.arange(locs.shape[0]):
        #a1 = numpy.random.uniform(-2.0,2.0)
        #a2 = numpy.random.uniform(-2.0,2.0)
        beamr = norm.pdf(x, loc=0, scale=numpy.random.uniform(lb,ub))#*norm.cdf(a1*x)
        beamp = norm.pdf(x, loc=0, scale=numpy.random.uniform(lb,ub))#*norm.cdf(a2*x)
        #beams_matrix[i,:] = beamr * numpy.exp(1j * beamp)
        beams_matrix[i,:] = beamr + 1j * beamp
        #plt.plot(x, beams_matrix[i,:].real)
        #plt.show()
        #plt.plot(x, beams_matrix[i,:].imag)
        #plt.show()
    
    return x,beams_matrix

# For simulating Electric Fields
def simulate_electric_fields_different_beams(locations_uvw,points,beam_matrices):
    elec_ss = numpy.zeros(shape=(locations_uvw.shape[0]),dtype=numpy.complex128)
    
    phase_source = numpy.zeros(shape=(points.shape[0]),dtype=numpy.complex128)
    print(phase_source.shape)
    for i in numpy.arange(points.shape[0]):
        phase_source[i] = numpy.exp(-1j*numpy.random.uniform(0,2*numpy.pi))
    
    for i in numpy.arange(locations_uvw.shape[0]):
        
        antenna_noise = numpy.random.normal(scale=4000.0)# + 1j * numpy.random.normal(scale=2000.0)
        for j,point in enumerate(points):
            l = point[0]
            m = point[1]
            n = numpy.sqrt(1.0-l**2-m**2) - 1.0
            bap = find_beam_at_point(beam_matrices[i,:],point)
            #print(bap)
            amp = point[2]
            #print(amp)
            u = locations_uvw[i,0]
            v = locations_uvw[i,1]
            w = locations_uvw[i,2]
            elec_ss[i] += amp*numpy.exp(-2j*numpy.pi*(u*l+v*m+w*n)) * bap * phase_source[j]
        #elec_ss[i] += antenna_noise
    return elec_ss

