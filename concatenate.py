import pyrap.tables
import lofar.parmdb
import os


msnames = ["L103929/L103929_SAP001_BAND0.MS", "L103931/L103931_SAP001_BAND0.MS", "L103933/L103933_SAP001_BAND0.MS", "L103935/L103935_SAP001_BAND0.MS", "L103937/L103937_SAP001_BAND0.MS", "L103939/L103939_SAP001_BAND0.MS"]

pyrap.tables.msutil.msconcat(msnames, "concatenated.MS")

pdb_concat_name = "concatenated.MS/ionosphere"
os.system("rm %s -rf" % pdb_concat_name)
pdb_concat = lofar.parmdb.parmdb(pdb_concat_name, create=True)

for msname in msnames:
    print msname
    pdb = lofar.parmdb.parmdb(msname + "/ionosphere")
    for parmname in pdb.getNames():
        print parmname
        v = pdb.getValuesGrid(parmname)
        pdb_concat.addValues(v)