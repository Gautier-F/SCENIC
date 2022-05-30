#!/bin/bash/

echo "Nom du dossier?"
read nom

mkdir -p "../${nom}"/{loomfiles,res,10x,scripts}
mkdir -p "../${nom}"/res/{fig,dat}

cp Readme.txt "../${nom}"/
cp \#2_SCN_grn.sh "../${nom}"/scripts
cp \#3_SCN_reg_pred.sh "../${nom}"/scripts
cp \#4_SCN_AUcells.sh "../${nom}"/scripts
cp script_scanpy.py "../${nom}"/scripts
cp \#1_qsub_scanpy.sh "../${nom}"/scripts
cp \#5_qsub_bin_mtx.sh "../${nom}"/scripts
cp \#6_qsub_Loom_to_Seurat.sh "../${nom}"/scripts
cp script_bin_mtx.py "../${nom}"/scripts
cp script_LoomToSeurat.R "../${nom}"/scripts
cp binarize.py "../${nom}"/scripts
cp ScSc.yaml "../${nom}"/scripts
cp config.yaml "../${nom}"/scripts
cp SnakeFile "../${nom}"/scripts
cp script_scanpy_snmk.py "../${nom}"/scripts