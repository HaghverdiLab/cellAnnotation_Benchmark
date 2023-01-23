#!/bin/bash
# -cwd

path=/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation
path_scripts="$path/Scripts"
path_data="$path/Data/Subsets/"
output="$path/Data/Predictions_curated/"
features="/fast/AG_Haghverdi/Carla_Moelbert/PaCMAP/Bonemarrow_genes.csv"
set=PBMC
test=$path/Data/Fulldata/$set"_Query"

for f in $features 1000 444; do
for set in PBMC; do #  Lung MotorCortex PBMC Kidney
  echo $set
  echo $test
  echo $f
  Filelist_length=$(ls -d $path_data/$set"10x_3090"* | wc -l)
  echo ${Filelist_length} 
  
  qsub -cwd -N Seurat  -l m_mem_free=32G -l h_rt=0:45:00 -t 1-${Filelist_length} \
  $path_scripts/Function_calls/seurat.sh $set $path_data $test $output/Seurat $f
  
  #qsub -cwd -N CellID  -l m_mem_free=64G -l h_rt=0:45:00  -t 1-${Filelist_length} \
  #$path_scripts/Function_calls/cellid.sh $set $path_data $test $output/CellID  $f
 
  #qsub -cwd -N SCN     -l m_mem_free=64G -l h_rt=0:45:00  -t 1-${Filelist_length} \
  #$path_scripts/Function_calls/scn.sh $set $path_data $test $output/SingleCellNet $f
  
  #qsub -cwd -N singleR -l m_mem_free=64G -l h_rt=0:45:00  -t 1-${Filelist_length} \
  #$path_scripts/Function_calls/singler.sh $set $path_data $test $output/SingleR $f

done
done