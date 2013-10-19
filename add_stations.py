from read_solutions import *
from fit_tec_per_source_pair import *
from calculate_piercepoints import *
import pyrap.measures
import utils


#No good solutions for RS409 RS310
bad_stations = ['RS409LBA', 'RS310LBA', 'RS106LBA', 'RS306LBA']
bad_stations = ['RS409LBA', 'RS310LBA']

beta = 5.0/3
datadir = "/home/vdtol/msss_peeling/L103929"
height = 200e3
order = 15

if 'source_names' not in locals():
    phases0,phases1, flags, mask, station_names, station_positions, source_names, source_positions, freqs, freqwidths, times, timewidths, pointing = read_solutions(datadir)
    station_position0 = station_positions[0,:]
    station_selection = find(sqrt(sum((station_positions - station_position0)**2, axis=1)) < 2e3)
    r0, source_selection,sols0,source_pairs = fit_tec_per_source_pair(phases0[:,station_selection,:,:], flags[:,station_selection,:,:], mask, freqs, propagate = True)
    r1, source_selection,sols1,source_pairs = fit_tec_per_source_pair(phases1[:,station_selection,:,:], flags[:,station_selection,:,:], mask, freqs, propagate = True)
    r = (r0+r1)/2


N_sources_selected = len(source_selection)
N_stations_selected = len(station_selection)
N_piercepoints = N_sources_selected * N_stations_selected
N_times = len(times)
N_stations = len(station_names)
N_sources = len(source_names)

D = resize( station_positions, ( N_stations, N_stations, 3 ) )
D = transpose( D, ( 1, 0, 2 ) ) - D
D = sqrt(sum( D**2, axis=2 ))

station_selection1 = station_selection


stations_to_add = array([i for i in range(len(station_names)) if i not in station_selection1 and station_names[i] not in bad_stations])
print station_names[stations_to_add]

q = r

while len(stations_to_add)>0:
    D1 = D[stations_to_add[:,newaxis], station_selection1[newaxis,:]]

    minimum_distance = amin(D1, axis=1)
    station_to_add = stations_to_add[argmin(minimum_distance)]
    print station_names[station_to_add]
    station_selection1 = append(station_selection1, station_to_add)
    N_stations_selected1 = len(station_selection1)
    
    # Remove station from list
    stations_to_add = stations_to_add[stations_to_add != station_to_add]
    
    sols_list = []
    eq_list = []
    min_e_list = []
    
    for ii, i in enumerate(source_selection):
        for jj, j in enumerate(source_selection):
            if j == i:
                break
            subband_selection = find(mask[i,:] * mask[j,:])
            if len(subband_selection)<=1:
                continue
            print '===================='
            print i,j, 'number of subbands = ', len(subband_selection)
            print '===================='
            p0 = phases0[i, station_selection1[:,newaxis], subband_selection[newaxis,:], :] - phases0[j, station_selection1[:,newaxis], subband_selection[newaxis,:], :]
            p0 = p0 - mean(p0, axis=0)[newaxis,:,:]
            p1 = phases1[i, station_selection1[:,newaxis], subband_selection[newaxis,:], :] - phases1[j, station_selection1[:,newaxis], subband_selection[newaxis,:], :]
            p1 = p1 - mean(p1, axis=0)[newaxis,:,:]
            A = zeros((len(subband_selection), 1))
            A[:,0] = 8.44797245e9/freqs[subband_selection]

            flags_source_pair = flags[i, station_selection1[:,newaxis], subband_selection[newaxis,:], :] * flags[j, station_selection1[:,newaxis], subband_selection[newaxis,:], :]
            constant_parms = zeros((1, N_stations_selected1), dtype = bool)
            sols = zeros((N_times, N_stations_selected1), dtype = float)
            
            for t_idx in range(N_times):
                min_e = Inf
                for offset in linspace(-0.1, 0.1,21) :
                    p_0 = zeros((1, N_stations_selected1), double)
                    p_0[0,:N_stations_selected1-1] = (q[ii, t_idx, :] - q[jj, t_idx, :])[newaxis,:]
                    p_0[0,-1] = offset
                
                    x = p0[:,:,t_idx].copy()
                    f = flags_source_pair[:,:,t_idx].copy()
                    sol0 = lofar.expion.baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                    sol0 -= mean(sol0)
                    residual = mod(dot(A,sol0) - x.T + pi, 2*pi) - pi
                    residual = residual[f.T==0]
                    e = var(residual)

                    x = p1[:,:,t_idx].copy()
                    f = flags_source_pair[:,:,t_idx].copy()
                    sol1 = lofar.expion.baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                    sol1 -= mean(sol1)
                    residual = mod(dot(A,sol1) - x.T + pi, 2*pi) - pi
                    residual = residual[f.T==0]
                    e += var(residual)

                    
                    if e<min_e:
                        min_e = e
                        sols[t_idx, :] = (sol0[0,:] + sol1[0,:])/2

                    #residual = mod(dot(A,sol) - x.T + pi, 2*pi) - pi
                    #residual = residual[f.T==0]
                    #ee = var(residual)
                    #e += ee
                    #e_list.append(sum(abs(residual) > 3*std(residual)))
                #print offset, e

            #weight = len(subband_selection) / min_e
            
            ### Remove outliers
            for kk in range(10):
                s = sols[:,-1].copy()
                selection = zeros(len(s), numpy.bool)
                for t_idx in range(len(s)):
                    start_idx = max(t_idx-10, 0)
                    end_idx = min(t_idx+10, len(s)) 
                    selection[t_idx] = sum(abs(s[start_idx:end_idx] - s[t_idx])<0.02) > (end_idx - start_idx - 8)
                outliers = find(logical_not(selection))
                if len(outliers) == 0:
                    break
                for t_idx in outliers:
                    print "outlier", t_idx
                    try:
                        idx0 = find(selection[:t_idx])[-1]
                    except IndexError:
                        idx0 = -1
                    try:
                        idx1 = find(selection[t_idx+1:])[0]+t_idx+1
                    except IndexError:
                        idx1 = -1
                    if idx0 == -1:
                        s[t_idx] = s[idx1]
                    elif idx1 == -1:
                        s[t_idx] = s[idx0]
                    else:
                        s[t_idx] = (s[idx0] * (idx1-t_idx) + s[idx1] * (t_idx-idx0)) / (idx1-idx0)
                
                    p_0 = zeros((1, N_stations_selected1), double)
                    p_0[0,:] = sols[t_idx,:]
                    p_0[0,-1] = s[t_idx]
                
                    x = p0[:,:,t_idx].copy()
                    f = flags_source_pair[:,:,t_idx].copy()
                    sol0 = lofar.expion.baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                    sol0 -= mean(sol0)

                    x = p1[:,:,t_idx].copy()
                    sol1 = lofar.expion.baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                    sol1 -= mean(sol1)
                    
                    print sols[t_idx, -1], (sol0[0,-1] + sol1[0,-1])/2, s[t_idx]
                    sols[t_idx, :] = (sol0[0,:] + sol1[0,:])/2
                    #sols[t_idx, -1] = s[t_idx]

            weight = 1.0
            sols_list.append(weight*sols)
            min_e_list.append(min_e)
            eq = zeros(N_sources)
            eq[ii] = weight
            eq[jj] = -weight
            eq_list.append(eq)
            
    sols = array(sols_list)
    B = array(eq_list)
    pinvB = pinv(B)

    q = dot(pinvB, sols.transpose([1,0,2]))
    #break
    if N_stations_selected1 == 28:
        break