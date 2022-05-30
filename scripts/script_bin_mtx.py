#### a utiliser sous env Scenic

import os 
import pandas as np
import loompy as lp
from binarize import *

con= lp.connect(snakemake.input[0], mode="r+", validate=False)

auc_mtx=pd.DataFrame(con.ca.RegulonsAUC, index = con.ca.CellID)

bin_mtx = binarize(auc_mtx)

mtx = np.array(bin_mtx[0])

con.ca.binary_mtx = mtx

con.close()
