#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -M fabien.gautier@univ-nantes.fr
#$ -m eba
pyscenic grn ../loomfiles/dat_filt_scenic.loom ../../Ressources/allTFs_hg38.txt -o ../res/dat/adj.csv --num_workers 20
