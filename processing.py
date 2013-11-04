#!/usr/bin/env python

import os
import numpy
from pylab import *
import lofar.parmdb

from read_solutions import *
from fit_tec_per_source_pair import *
from add_stations import add_stations
from calculate_piercepoints import *
from fit_screen_to_tec import *
from make_images import *
#from plot_images import *

for obsid in ["L103929", "L103931", "L103933", "L103935", "L103937", "L103939"] :
    datadir = "/data/scratch/vdtol/MSSS_LBA/data/%s" % obsid
    height = 200e3
    order = 15

    iondatafile = 'iondata-%s.npz' % obsid
    
    if os.path.exists(iondatafile):
        iondata = numpy.load( iondatafile )

        phases0 = iondata['phases0']
        phases1 = iondata['phases1']
        flags = iondata['flags']
        mask = iondata['mask']
        station_names = iondata['station_names']
        station_positions = iondata['station_positions']
        station_selection = iondata['station_selection']
        source_names = iondata['source_names']
        source_positions = iondata['source_positions']
        source_selection = iondata['source_selection']
        freqs = iondata['freqs']
        freqwidths = iondata['freqwidths']
        times = iondata['times']
        timewidths = iondata['timewidths']
        pointing = iondata['pointing']
        r = iondata['r']
        sols = iondata['sols']

    else:
        # read the direction dependent gain solutions from the parmdbs
        phases0,phases1, flags, mask, station_names, station_positions, source_names, source_positions, freqs, freqwidths, times, timewidths, pointing = read_solutions(datadir)

        # select the station that are less then 2km from the first station (CS001)
        station_position0 = station_positions[0,:]
        station_selection = find(sqrt(sum((station_positions - station_position0)**2, axis=1)) < 2e3)

        # fit a tec value to the phase solutions 
        r0, source_selection,sols0,source_pairs = fit_tec_per_source_pair(phases0[:,station_selection,:,:], flags[:,station_selection,:,:], mask, freqs, propagate = True)
        r1, source_selection,sols1,source_pairs = fit_tec_per_source_pair(phases1[:,station_selection,:,:], flags[:,station_selection,:,:], mask, freqs, propagate = True)

        # take the mean of the two polarizations
        r = (r0+r1)/2
        
        station_selection, sols, r = add_stations(station_selection, phases0,phases1, flags, mask, station_names, station_positions, source_names, source_selection, times,freqs, r)    

        numpy.savez( 'iondata-%s.npz' % obsid, 
            phases0 = phases0,
            phases1 = phases1,
            flags = flags,
            mask = mask,
            station_names = station_names,
            station_positions = station_positions,
            source_names = source_names,
            source_positions = source_positions,
            freqs = freqs,
            freqwidths = freqwidths,
            times = times,
            timewidths = timewidths,
            pointing = pointing,
            sols = sols,
            r = r,
            station_selection = station_selection,
            source_pairs = source_pairs,
            source_selection = source_selection)

    
    r = r[source_selection, :, :]
    
    N_sources = len(source_selection)
    N_times = len(times)
    N_stations = len(station_selection)

    N_piercepoints = N_sources * N_stations
    
    print N_sources, N_times, N_stations, N_piercepoints, r.shape
    
    rr = reshape(r.transpose([0,2,1]), [ N_piercepoints, N_times])


    pp = calculate_piercepoints(station_positions[station_selection,:], source_positions[source_selection,:], times, height)

    # fit a TEC screen

    for subbandnr in [0]:
        msname = datadir + "/%s_SAP001_BAND%i.MS" % (obsid, subbandnr)
        imagedir = "/data/scratch/vdtol/MSSS_LBA/images/%s/BAND%i" % (obsid, subbandnr)

        parms, tec_fit = fit_screen_to_tec(station_names[station_selection,:], source_names[source_selection], pp, rr, times, timewidths, height, order, freqs[[subbandnr]], freqwidths[[subbandnr]])

        # write the TEC screen to a parmdb
        os.system('rm %s/ionosphere/ -rf' % msname)
        pdb = lofar.parmdb.parmdb("%s/ionosphere/" % msname, create=True)
        pdb.addValues(parms)
        pdb.addDefValues({'Gain:0:0:Ampl': pdb.makeDefValue(1.0), 'Gain:1:1:Ampl': pdb.makeDefValue(1.0)})
        pdb.flush()

        os.system("calibrate-stand-alone --no-columns --parmdb-name ionosphere  %s bbs.parset none" % msname)

        try:
            os.makedirs(imagedir)
        except OSError:
            pass


        make_images(msname, imagedir, source_names)

        #plot_images(imagedir, source_names, source_positions)
