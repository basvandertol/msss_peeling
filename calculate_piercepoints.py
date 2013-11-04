import pyrap.measures
import utils
import numpy

def calculate_piercepoints(station_positions, source_positions, times, height = 200e3):

    N_sources = source_positions.shape[0]
    N_stations = station_positions.shape[0]
    N_piercepoints = N_stations * N_sources
    N_times = len(times)
    
    me = pyrap.measures.measures()
    position = me.position('ITRF', '%fm' % station_positions[0,0], '%fm' % station_positions[0,1], '%fm' % station_positions[0,2])
    me.doframe(position)


    pp = numpy.zeros((N_times, N_piercepoints,3))

    for k in range(N_times):
        epoch = me.epoch('UTC', '%fs' % times[k])
        me.doframe(epoch)
        pp_idx = 0
        for i in range(N_sources):
            ra = source_positions[i,0]
            dec = source_positions[i,1]
            d = me.direction('J2000', '%frad' % ra, '%frad' % dec)
            d1 = me.measure(d, 'ITRF')
            phi = d1['m0']['value']
            theta = d1['m1']['value']
            dx = numpy.cos(theta)*numpy.cos(phi)
            dy = numpy.cos(theta)*numpy.sin(phi)
            dz = numpy.sin(theta)
            direction = numpy.array((dx,dy,dz))
            for station_position in station_positions:
                pp[k, pp_idx, :] = utils.calc_piercepoint(station_position, direction, height)
                pp_idx += 1

    
    return pp
