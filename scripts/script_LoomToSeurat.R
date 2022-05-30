#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)

# test if arguments are indicated: if not, return an error
#if (length(args)==0) {
# stop("project name, path to loomfile, and path to 10x must be supplied", call.=FALSE)
#}

library(loomR)
library(Seurat)

project_name= list.files("../10x/") #args[1]
path_to_loom = "../res/dat/bin_mtx_SCENIC.loom" #args[2]
path_to_10x = paste("../10x/", project_name, sep="") #args[3]


con = connect(path_to_loom, skip.validate= TRUE)

cell_id = con$col.attrs$CellID
C_ID = cell_id[1:cell_id$dims]

regAUC = con$col.attrs$RegulonsAUC
df_regAUC = regAUC[1]
for(i in 2:regAUC$dims){
  df_regAUC[i,]= regAUC[i][1,]
}

Bin_mtx= con$col.attrs$binary_mtx
df_bin=Bin_mtx[1:Bin_mtx$dims[1],1:Bin_mtx$dims[2]] 
df_bin = as.data.frame(t(df_bin))
colnames(df_bin)=colnames(df_regAUC)
rownames(df_bin)=C_ID


#inclusion dans Seurat object
Seurat=Read10X(path_to_10x)
Seurat=CreateSeuratObject(counts=Seurat, project=project_name, min.cells = 3,
                          min.features=200)


Seurat@meta.data$in_SCN = colnames(Seurat)%in%C_ID
Seurat <- subset(Seurat, subset = in_SCN)

Seurat@meta.data$RegAUC=df_regAUC
Seurat@meta.data$Bin_mtx=df_bin

Regulons = con$row.attrs$Regulons
Reg = Regulons$read(args= list(1: Regulons$dims))
genes =  con$row.attrs$Gene
Genes = genes$read(args= list(1:genes$dims[1]))
rownames(Reg) = Genes
Seurat[["RNA"]][["Regulons"]]=Reg
con$close_all()

saveRDS(Seurat, paste("../res/dat/", project_name, ".rds", sep = ""))





