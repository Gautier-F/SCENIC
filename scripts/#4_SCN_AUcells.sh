#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -M fabien.gautier@univ-nantes.fr
#$ -m eba

pyscenic aucell \
../loomfiles/dat_filt_scenic.loom \
../res/dat/reg.csv \
--output ../res/dat/bin_mtx_SCENIC.loom
--num_workers 20 

