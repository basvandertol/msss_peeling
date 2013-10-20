import utils
from pylab import *
import lofar.expion


def fit_tec_per_source_pair(phases, flags, mask, freqs, init_sols = None, init_sols_per_pair = False, propagate = False):

    sols_list = []
    eq_list = []

    N_sources = phases.shape[0]
    N_stations = phases.shape[1]
    N_times = phases.shape[3]
    
    if init_sols is None:
        init_sols = zeros((N_sources, N_times, N_stations), dtype = float)

    source_pairs = []
    
    k = 0
    for i in range(N_sources):
        for j in range(i):
            subband_selection = find(mask[i,:] * mask[j,:])
            if len(subband_selection)<=1:
                continue
            
             
            print i,j
            source_pairs.append((i,j))
            p = phases[i, :, subband_selection, :] - phases[j, :, subband_selection, :]
            A = zeros((len(subband_selection), 1))
            A[:,0] = 8.44797245e9/freqs[subband_selection]

            flags_source_pair = flags[i, :, subband_selection, :] * flags[j, :, subband_selection, :]
            constant_parms = zeros((1, N_stations), dtype = bool)
            sols = zeros((N_times, N_stations), dtype = float)
            
            if init_sols_per_pair:
                p_0 = init_sols[k, 0, :][newaxis,:]
            else:
                p_0 = (init_sols[i, 0, :] - init_sols[j, 0, :])[newaxis,:]
            
            for t_idx in range(N_times):
                x = p[:,:,t_idx].copy()
                f = flags_source_pair[:,:,t_idx].copy()
                if not propagate:
                    if init_sols_per_pair:
                        p_0 = init_sols[k, t_idx, :][newaxis,:]
                    else:
                        p_0 = (init_sols[i, t_idx, :] - init_sols[j, 0, :])[newaxis,:]
                sol = lofar.expion.baselinefitting.fit(x, A, p_0, f, constant_parms)
                if propagate:
                    p_0 = sol.copy()
                sols[t_idx, :] = sol[0,:]
            #subplot(1,2,2)
            sols = sols[:, :]-mean(sols[:, :], axis=1)[:,newaxis]
            weight = len(subband_selection)
            sols_list.append(sols*weight)
            eq = zeros(N_sources)
            eq[i] = weight
            eq[j] = -weight
            eq_list.append(eq)
            k += 1
                
    sols = array(sols_list)
    B = array(eq_list)
    source_selection = find(sum(abs(B), axis=0))
    N_sources = len(source_selection)


    pinvB = pinv(B)

    r = dot(pinvB, sols.transpose([1,0,2]))

    r = r[source_selection,:,:]

    return r, source_selection, sols, source_pairs
    
    #r = r - r[:,:,[0]]
    #r = r - r[[0],:,:]
