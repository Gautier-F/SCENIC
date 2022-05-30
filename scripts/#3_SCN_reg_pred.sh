#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -M fabien.gautier@univ-nantes.fr
#$ -m eba

pyscenic ctx ../res/dat/adj.csv \
../../Ressources/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname ../../Ressources/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname ../loomfiles/dat_filt_scenic.loom \
--mode "dask_multiprocessing" \
--output ../res/dat/reg.csv \
--num_workers 20 \
--mask_dropouts 
