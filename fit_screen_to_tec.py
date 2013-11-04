import pyrap.measures
import utils
import os
from pylab import *

def fit_screen_to_tec(station_names, source_names, pp, rr, times, timewidths, height, order, freqs0, freqwidths0):

    N_stations = len(station_names)
    N_sources = len(source_names)
    N_times = len(times)
    
    tec_fit_all = zeros((N_times, N_sources, N_stations))

    A = concatenate([kron(eye(N_sources), ones((N_stations,1))), kron(ones((N_sources,1)), eye(N_stations))], axis=1)

    N_piercepoints = N_sources*N_stations
    P = eye(N_piercepoints) - dot(dot(A, pinv(dot(A.T, A))), A.T)

    freqs = freqs0
    freqwidths = freqwidths0
    
    N_freqs = 1

    parms = {}    
    v = {}
    v['times'] = times
    v['timewidths'] = timewidths 
    v['freqs'] = freqs
    v['freqwidths'] = freqwidths

    for station_name in station_names:
        
        #v['values'] = zeros((N_times, N_freqs), dtype=double)
        #parmname = 'Gain:0:0:Phase:' + station_name
        #parms[parmname] = v.copy()
        
        #v['values'] = ones((N_times, N_freqs), dtype=double)
        #parmname = 'Gain:0:0:Ampl:' + station_name
        #parms[parmname] = v.copy()
        
        #v['values'] = zeros((N_times, N_freqs), dtype=double)
        #parmname = 'Gain:1:1:Phase:' + station_name
        #parms[parmname] = v.copy()
        
        #v['values'] = ones((N_times, N_freqs), dtype=double)
        #parmname = 'Gain:1:1:Ampl:' + station_name
        #parms[parmname] = v.copy()
        
        for source_name in source_names:
            
            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'Piercepoint:X:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()
            
            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'Piercepoint:Y:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()
            
            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'Piercepoint:Z:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()
            
            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'TECfit_white:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()
            
            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'TECfit_white:0:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()
            
            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'TECfit_white:1:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()
            
    for k in range(N_times):
        D = resize( pp[k,:,:], ( N_piercepoints, N_piercepoints, 3 ) )
        D = transpose( D, ( 1, 0, 2 ) ) - D
        D2 = sum( D**2, axis=2 )
        beta = 5.0/3
        r_0 = 10e3
        C = -(D2 / ( r_0**2 ) )**( beta / 2.0 )/2.0
        P1 = eye(N_piercepoints) - ones((N_piercepoints, N_piercepoints)) / N_piercepoints
        C1 = dot(dot(P1, C ), P1)
        U,S,V = svd(C1)

        B = dot(P, U[:,:order])
        pinvB = pinv(B,rcond=1e-3)

        rr1 = dot(P, rr[:,k])

        tec_fit = dot(U[:,:order],dot(pinvB, rr1))
        tec_fit_all[k,:,:] = tec_fit.reshape((N_sources, N_stations))

        tec_fit_white = dot(pinv(C), tec_fit)
        pp_idx = 0
        for source_name in source_names:
            for station_name in station_names:
                
                parmname = 'Piercepoint:X:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k,0] = pp[k,pp_idx,0]
                
                parmname = 'Piercepoint:Y:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k,0] = pp[k,pp_idx,1]
                
                parmname = 'Piercepoint:Z:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k,0] = pp[k,pp_idx,2]
                
                parmname = 'TECfit_white:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k,0] = tec_fit_white[pp_idx]
                
                parmname = 'TECfit_white:0:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k,0] = tec_fit_white[pp_idx]
                
                parmname = 'TECfit_white:1:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k,0] = tec_fit_white[pp_idx]
                
                pp_idx += 1

                
    time_start = times[0] - timewidths[0]/2
    time_end = times[-1] + timewidths[-1]/2

    v[ 'times' ] = array([(time_start + time_end) / 2])
    v[ 'timewidths' ] = array([time_end - time_start])
                
    v_r0 = v.copy()
    v_r0['values'] = array(r_0, dtype=double, ndmin=2)
    parms['r_0'] = v_r0

    v_beta = v.copy()
    v_beta['values'] = array(beta, dtype=double, ndmin=2)
    parms['beta'] = v_beta

    v_height = v.copy()
    v_height['values'] = array(height, dtype=double, ndmin=2)
    parms['height'] = v_height
    
    return parms, tec_fit_all