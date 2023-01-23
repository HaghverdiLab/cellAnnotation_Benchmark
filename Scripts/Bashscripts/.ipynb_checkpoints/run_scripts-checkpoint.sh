#!/bin/bash
# -cwd

path=/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation
path_scripts="$path/Scripts"
path_data="$path/Data/Subsets_runtime/"
output="$path/Data/Predictions/"


for set in PBMC; do 
  echo $set
  test=$path/Data/Fulldata/$set"_Query"
  echo $test
  
  Filelist_length=$(ls -d $path_data/$set* | wc -l)
  echo ${Filelist_length} 
  qsub -cwd -N ItClust  -l m_mem_free=16G -l h_rt=0:15:00 -l cpu=8  -t 1-${Filelist_length}  \
  $path_scripts/Function_calls/itClust.sh $set $path_data $test $output/ItClust 
  
  qsub -cwd -N Seurat  -l m_mem_free=16G -l h_rt=0:15:00  -l cpu=8 -t 1-${Filelist_length} \
  $path_scripts/Function_calls/seurat.sh $set $path_data $test $output/Seurat
  
  qsub -cwd -N CellID  -l m_mem_free=16G -l h_rt=0:15:00  -l cpu=8  -t 1-${Filelist_length} \
  $path_scripts/Function_calls/cellid.sh $set $path_data $test $output/CellID 
 
  qsub -cwd -N SCN   -l m_mem_free=16G -l h_rt=0:15:00 -l cpu=8  -t 1-${Filelist_length} \
  $path_scripts/Function_calls/scn.sh $set $path_data $test $output/SingleCellNet
  
  qsub -cwd -N singleR -l m_mem_free=16G -l h_rt=0:15:00 -l cpu=8 -t 1-${Filelist_length} \
  $path_scripts/Function_calls/singler.sh $set $path_data $test $output/SingleR

done