import xlsxwriter
import pandas as pd
import loompy as lp
import sys


loom  = sys.argv[1]
file = sys.argv[2]

con = lp.connect(loom, mode = "r+", validate= False)
bin_mtx=pd.DataFrame(con.ca.binary_mtx, index= con.ca.CellID, columns=pd.DataFrame(con.ca.RegulonsAUC).columns)
hdr = list(bin_mtx.columns)
con.close()

writer = pd.ExcelWriter(file, engine="xlsxwriter")
bin_mtx.to_excel(writer, sheet_name='Bmtx', header= hdr)
writer.save()
