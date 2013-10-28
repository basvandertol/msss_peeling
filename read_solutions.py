#!/usr/bin/env python

import glob
import lofar.parmdb
from pylab import *
import pyrap.tables
import pyrap.measures
import lofar.expion.baselinefitting
import time
import os

import utils


def read_solutions(datadir = 'L103929'):
    mslist = sort(glob.glob(datadir + '/*.MS'))
    N_freqs = len(mslist)
    N_times = 66

    t_ant = pyrap.tables.table(mslist[0] + '/ANTENNA')
    N_ant = len(t_ant)

    s = set() 
    for i, msname in enumerate(mslist[:]):
        instrumentname = msname + '/instrument'
        #print instrumentname
        pdb = lofar.parmdb.parmdb( instrumentname )
        s = s.union(utils.get_source_list( pdb, ["*"] ))
    source_names = array(sorted(s))
    N_sources = len(source_names)


    m = zeros((N_sources, N_freqs))
    for i, msname in enumerate(mslist[:]):
        instrumentname = msname + '/instrument'
        #print instrumentname
        pdb = lofar.parmdb.parmdb( instrumentname )
        source_list1 = utils.get_source_list( pdb, ["*"] )
        for source in source_list1:
            if source in source_names:
                m[find(source_names == source), i] = 1
            
    skydb = lofar.parmdb.parmdb(mslist[0] + '/sky')


    instrumentname = mslist[0] + '/instrument'
    pdb = lofar.parmdb.parmdb( instrumentname )
    parmname = 'DirectionalGain:0:0:Real:CS001LBA' + ':' + source_names[2]
    times = pdb.getValuesGrid(parmname)[parmname]['times']
    timewidths = pdb.getValuesGrid(parmname)[parmname]['timewidths']

    station_positions = t_ant.getcol("POSITION")
    station_names = array(t_ant.getcol("NAME"))
    N_stations = len(station_names)

    phases = zeros((N_sources, N_stations, N_freqs, N_times))
    phases1 = zeros((N_sources, N_stations, N_freqs, N_times))
    flags = ones((N_sources, N_stations, N_freqs, N_times))

    freqs = zeros(N_freqs)
    freqwidths = zeros(N_freqs)
    
    source_positions = zeros((N_sources,2))

    for i, source1 in enumerate(source_names):

        parmname = "Dec:" + source1
        dec1 = skydb.getDefValues(parmname)[parmname][0,0]
        parmname = "Ra:" + source1
        ra1 = skydb.getDefValues(parmname)[parmname][0,0]
        source_positions[i, ] = [ra1, dec1]
        
        for k, msname in enumerate(mslist[:]):
            if m[i, k]:
                instrumentname = msname + '/instrument'
                pdb = lofar.parmdb.parmdb( instrumentname )
                parmname = 'DirectionalGain:0:0:Real:CS001LBA' + ':' + source1
                freqs[k] = pdb.getValuesGrid(parmname)[parmname]['freqs'][0]
                freqwidths[k] = pdb.getValuesGrid(parmname)[parmname]['freqwidths'][0]
                for l, station in enumerate(station_names):
                    parmname = 'DirectionalGain:0:0:Real:' + station + ':' + source1
                    v1_re = pdb.getValuesGrid(parmname)[parmname]['values']
                    
                    parmname = 'DirectionalGain:0:0:Imag:' + station + ':' + source1
                    v1_im = pdb.getValuesGrid(parmname)[parmname]['values']
                    v1 = arctan2( v1_im, v1_re)
                    phases[i, l, k,:] = v1[:,0]
                    flags[i, l, k,:] = (v1_re[:,0] == 1.0)

                    parmname = 'DirectionalGain:1:1:Real:' + station + ':' + source1
                    v1_re = pdb.getValuesGrid(parmname)[parmname]['values']
                    parmname = 'DirectionalGain:1:1:Imag:' + station + ':' + source1
                    v1_im = pdb.getValuesGrid(parmname)[parmname]['values']
                    v1 = arctan2( v1_im, v1_re)
                    phases1[i, l, k,:] = v1[:,0]

    mask = m
    
    t = pyrap.tables.table(mslist[0])
    t_field = pyrap.tables.table(t.getkeyword('FIELD'))
    ra = t_field[0]['PHASE_DIR'][0][0]
    dec = t_field[0]['PHASE_DIR'][0][1]
    pointing = array([cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)])


    return phases, phases1, flags, mask, station_names, station_positions, source_names, source_positions, freqs, freqwidths, times,timewidths, pointing 
