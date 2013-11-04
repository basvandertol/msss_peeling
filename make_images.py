#!/usr/bin/env python

import pyrap.images
import lofar.parameterset
import os

UVmax = 1.0
cellsize = '40arcsec'


def make_images(msname, imagedir, source_names) :
    
    #Make an image of the original data
    
    #ps = lofar.parameterset.parameterset('bbs.parset.template')
    #ps.replace('Strategy.Steps', "[correct]")
    #ps.writeFile('bbs.parset.tmp')
    #os.system('calibrate-stand-alone  --no-columns --parmdb-name ionosphere %s bbs.parset.tmp none' % msname)
    #os.system('rm %s/original.residual -rf' % imagedir)
    #os.system('rm %s/original.model -rf'  % imagedir )
    #os.system('awimager ms=%s data=CORRECTED_DATA image=%s/original operation=csclean niter=100 UVmax=%f cellsize=%s npix=2000' % (msname, imagedir, UVmax, cellsize))

    # Make image per source corrected for that source
    
    #for source_name in list(source_names)[1:]:
        #print source_name
        #ps = lofar.parameterset.parameterset('bbs.parset.template')
        #ps.replace('Step.correct2.Model.Sources', "[" + source_name + "]")
        #ps.writeFile('bbs.parset.tmp')
        #os.system('calibrate-stand-alone  --no-columns --parmdb-name ionosphere %s bbs.parset.tmp none' % msname)
        #os.system('rm %s/%s.residual -rf' % (imagedir, source_name))
        #os.system('rm %s/%s.model -rf' % (imagedir, source_name))
        #os.system('awimager ms=%s data=CORRECTED_DATA image=%s/%s-ampl operation=csclean niter=100 UVmax=%f cellsize=%s npix=2000' % (msname, imagedir, source_name, UVmax, cellsize))

    # make image with A-projection
    
    ps = lofar.parameterset.parameterset('bbs.parset.template')
    ps.replace('Step.correct2.Model.Sources', "[]")
    ps.writeFile('bbs.parset.tmp')
    os.system('calibrate-stand-alone --no-columns --parmdb-name ionosphere %s bbs.parset.tmp none' % msname)
    os.system('rm %s/aprojection.residual -rf' % imagedir)
    os.system('rm %s/aprojection.model -rf'  % imagedir)
    os.system('awimager ms=%s data=CORRECTED_DATA image=%s/aprojection operation=csclean niter=100 UVmax=%f cellsize=%s npix=2000 applyIonosphere=1 timewindow=10 SpheSupport=45 parmdbname=ionosphere threshold=20Jy' % (msname, imagedir, UVmax, cellsize))
        
