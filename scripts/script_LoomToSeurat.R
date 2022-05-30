
if (!require(loomR)){
	devtools::install_github(repo="mojaveazure/loomR",
							ref="develop")
}							

library(loomR)
library(Seurat)

project_name = Sys.getenv("DIR_NAME") 
path_to_loom = snakemake@input[["loom"]] 
path_to_10x = paste("../../../../matrix_10x/", project_name, sep = "") 


con = connect(path_to_loom, skip.validate=TRUE)

cell_id = con$col.attrs$CellID
C_ID = cell_id[1:cell_id$dims]

regAUC = con$col.attrs$RegulonsAUC
df_regAUC = regAUC[1]
for(i in 2:regAUC$dims){
  df_regAUC[i,] = regAUC[i][1,]
}

Bin_mtx = con$col.attrs$binary_mtx
df_bin = Bin_mtx[1:Bin_mtx$dims[1],1:Bin_mtx$dims[2]] 
df_bin = as.data.frame(t(df_bin))
colnames(df_bin) = colnames(df_regAUC)
rownames(df_bin) = C_ID


#inclusion dans Seurat object
Seurat = Read10X(path_to_10x)
Seurat = CreateSeuratObject(counts=Seurat, project=project_name, min.cells=3,
                          min.features=200)


Seurat@meta.data$in_SCN = colnames(Seurat)%in%C_ID
Seurat <- subset(Seurat, subset=in_SCN)

Seurat@meta.data$RegAUC = df_regAUC
Seurat@meta.data$Bin_mtx = df_bin

Regulons = con$row.attrs$Regulons
Reg = Regulons$read(args=list(1: Regulons$dims))
genes =  con$row.attrs$Gene
Genes = genes$read(args=list(1:genes$dims[1]))
rownames(Reg) = Genes
Seurat[["RNA"]][["Regulons"]] = Reg
con$close_all()

saveRDS(Seurat, paste("res/dat/", project_name, ".rds", sep=""))





