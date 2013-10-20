from pylab import *
import os


def calc_piercepoint(pos, direction, height):
   pp = zeros(3)
   earth_ellipsoid_a = 6378137.0;
   earth_ellipsoid_a2 = earth_ellipsoid_a*earth_ellipsoid_a;
   earth_ellipsoid_b = 6356752.3142;
   earth_ellipsoid_b2 = earth_ellipsoid_b*earth_ellipsoid_b;
   earth_ellipsoid_e2 = (earth_ellipsoid_a2 - earth_ellipsoid_b2) / earth_ellipsoid_a2;
 
   ion_ellipsoid_a = earth_ellipsoid_a + height;
   ion_ellipsoid_a2_inv = 1.0 / (ion_ellipsoid_a * ion_ellipsoid_a);
   ion_ellipsoid_b = earth_ellipsoid_b + height;
   ion_ellipsoid_b2_inv = 1.0 / (ion_ellipsoid_b * ion_ellipsoid_b);

   x = pos[0]/ion_ellipsoid_a;
   y = pos[1]/ion_ellipsoid_a;
   z = pos[2]/ion_ellipsoid_b;
   c = x*x + y*y + z*z - 1.0;
    
   dx = direction[0] / ion_ellipsoid_a
   dy = direction[1] / ion_ellipsoid_a
   dz = direction[2] / ion_ellipsoid_b
   a = dx*dx + dy*dy + dz*dz
   b = x*dx + y*dy  + z*dz
   alpha = (-b + sqrt(b*b - a*c))/a
   pp = pos[0] + alpha*direction
   normal_x = pp[0] * ion_ellipsoid_a2_inv
   normal_y = pp[1] * ion_ellipsoid_a2_inv
   normal_z = pp[2] * ion_ellipsoid_b2_inv
   norm_normal2 = normal_x*normal_x + normal_y*normal_y + normal_z*normal_z
   norm_normal = sqrt(norm_normal2)
   sin_lat2 = normal_z*normal_z / norm_normal2
   g = 1.0 - earth_ellipsoid_e2*sin_lat2
   sqrt_g = sqrt(g)
   M = earth_ellipsoid_b2 / ( earth_ellipsoid_a * g * sqrt_g )
   N = earth_ellipsoid_a / sqrt_g

   local_ion_ellipsoid_e2 = (M-N) / ((M+height)*sin_lat2 - N - height);
   local_ion_ellipsoid_a = (N+height) * sqrt(1.0 - local_ion_ellipsoid_e2*sin_lat2)
   local_ion_ellipsoid_b = local_ion_ellipsoid_a*sqrt(1.0 - local_ion_ellipsoid_e2)

   z_offset = ((1.0-earth_ellipsoid_e2)*N + height - (1.0-local_ion_ellipsoid_e2)*(N+height)) * sqrt(sin_lat2)

   x1 = pos[0]/local_ion_ellipsoid_a
   y1 = pos[1]/local_ion_ellipsoid_a
   z1 = (pos[2]-z_offset)/local_ion_ellipsoid_b
   c1 = x1*x1 + y1*y1 + z1*z1 - 1.0

   dx = direction[0] / local_ion_ellipsoid_a
   dy = direction[1] / local_ion_ellipsoid_a
   dz = direction[2] / local_ion_ellipsoid_b
   a = dx*dx + dy*dy + dz*dz
   b = x1*dx + y1*dy  + z1*dz
   alpha = (-b + sqrt(b*b - a*c1))/a

   pp = pos + alpha*direction
   return pp
   


def get_source_list( pdb, source_pattern_list ):
   source_list = []
   for pattern in source_pattern_list :
      parmname_list = pdb.getNames( 'DirectionalGain:?:?:*:*:' + pattern )
      source_list.extend([n.split(':')[-1] for n in parmname_list])
      parmname_list = pdb.getNames( 'RotationAngle:*:' + pattern )
      source_list.extend([n.split(':')[-1] for n in parmname_list])
   return sorted(set(source_list))

def get_station_list( pdb, station_pattern_list, DirectionalGainEnable ):
   station_list = []
   for pattern in station_pattern_list :
      parmname_list = pdb.getNames( { True : 'DirectionalGain:?:?:*:'+pattern+':*', False: 'Gain:?:?:*:' + pattern}[DirectionalGainEnable] )
      station_list.extend(sorted(set([n.split(':')[{True : -2, False : -1}[DirectionalGainEnable]] for n in parmname_list])))
   return station_list

def makemovie(sols, p, name = 'movie'):
    low = min(sols.flatten())
    high = max(sols.flatten())
    
    lower = amin(p, axis=0)
    upper = amax(p, axis=0)
    extent = upper - lower

    lower = lower - 0.05 * extent
    upper = upper + 0.05 * extent
    
    f = figure(figsize = (7,7))    
    for i in range(sols.shape[0]):
        clf()
        sm = cm.ScalarMappable(cmap = cm.jet,norm=normalize(vmin = low, vmax=high))
        sm._A = []
        gca().set_aspect('equal')
        title(str(i))
        colorbar(sm)    
        print sols[i,0]
        for j in range(sols.shape[1]-1,-1,-1):
            plot(p[j,0]*1e-3, p[j,1]*1e-3, {True: 's', False:'o'}[j==-1], markerfacecolor = sm.to_rgba(sols[i,j]), markersize=10)
        xlabel("x (km)")
        ylabel("y (km)")
        xlim(lower[0]*1e-3, upper[0]*1e-3)
        ylim(lower[1]*1e-3, upper[1]*1e-3)
        savefig('frame%0.3i.png' % i)
    #os.system("mencoder mf://frame???.png -o %s.mpeg -mf type=png:fps=3  -ovc x264 -ffourcc DX50 -noskip -oac copy" % name)
    os.system("ffmpeg -y -r 3 -i frame%03d.png  -c:v libvpx -b:v 1M -c:a libvorbis " + name + ".webm")
    close(f)
   

def plot_source_positions( source_positions, source_names, pointing, pairs = ())   :
    
    x,y,z = pointing
    east = array( [-y, x, 0])
    east = east / norm(east)

    north = array([ -x, -y, (x*x + y*y)/z])
    north = north/norm(north)

    T = concatenate([east[:,newaxis], north[:,newaxis]], axis=1)
    
    for i in range(len(source_names)):
        ra = source_positions[i,0]
        dec = source_positions[i,1]
        vv = array([cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)])
        q = dot(vv, T)
        plot(-q[0], q[1], 'o', markerfacecolor='k')
        text(-q[0], q[1], source_names[i])
        
    for pair in pairs:
        ra = source_positions[pair[0],0]
        dec = source_positions[pair[0],1]
        vv = array([cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)])
        q0 = dot(vv, T)
        ra = source_positions[pair[1],0]
        dec = source_positions[pair[1],1]
        vv = array([cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)])
        q1 = dot(vv, T)
        plot([-q0[0],-q1[0]], [q0[1],q1[1]], 'k')
        

def fit_gradient(station_positions, station_selection_core, source_names, sols, source_pairs, times):
    
    N_sources = len(source_names)
    
    station_position0 = station_positions[0,:]

    x,y,z = station_position0
    east = array( [-y, x, 0])
    east = east / norm(east)

    north = array([ -x, -y, (x*x + y*y)/z])
    north = north/norm(north)

    up = station_position0
    up = up/norm(up)

    T = concatenate([east[:,newaxis], north[:,newaxis]], axis=1)
    p = station_positions
    p1 = dot(p[:,:], T)

    s = sols.transpose((0,2,1))
    print s[:,station_selection_core,:].shape
    print pinv(p1[station_selection_core,:]).shape
    print dot(pinv(p1[station_selection_core,:]), s[:,station_selection_core,:]).shape
    
    residual = s - dot(p1,dot(pinv(p1[station_selection_core,:]), s[:,station_selection_core,:]).transpose((1,0,2))).transpose((1,0,2))
    residual = residual.transpose((0,2,1))
    #imshow(abs(residual), interpolation='nearest')
    #hist(residual.flatten(),100)
    
    return residual