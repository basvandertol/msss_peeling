import pyrap.tables
import os


msnames = ["L103929/L103929_SAP001_BAND0.MS", "L103931/L103931_SAP001_BAND0.MS", "L103933/L103933_SAP001_BAND0.MS", "L103935/L103935_SAP001_BAND0.MS", "L103937/L103937_SAP001_BAND0.MS", "L103939/L103939_SAP001_BAND0.MS"]


threshold = 600

for msname in msnames:
    t = pyrap.tables.table(msname, readonly = False)
    if os.path.exists(msname + '.flags'):
        t_flags = pyrap.tables.table(msname + '.flags')
    else:
        t_flags = t.select("FLAG").copy(msname + '.flags')
        
    f = t_flags.getcol("FLAG")
    d = t.getcol("CORRECTED_DATA")
    
    f = logical_or(f, abs(d)>threshold)
    t.putcol("FLAG", f)
    t.flush()
    