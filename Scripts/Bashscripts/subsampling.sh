#!/bin/sh
source /home/cmoelbe/bin/anaconda3/bin/activate r-environment

path=/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation
path_scripts=$path/Scripts
output="$path/Data/Subsets_runtime"
data=$path/Data/Fulldata/


#qsub -N balanced -l m_mem_free=32G -pe smp 2 -l h_rt=08:00:00 \
#$path_scripts/Function_calls/subsampling.sh $data/PBMC10x_Reference $output 300 \
#250,500,1000 $path_scripts 
    
#qsub -N bootstrap -l m_mem_free=32G -pe smp 2 -l h_rt=08:00:00 \
#$path_scripts/Function_calls/subsampling.sh $data/PBMC10x_Reference $output 300 \
#100 $path_scripts    

#qsub -N mosaic -l m_mem_free=32G -pe smp 2 -l h_rt=08:00:00 \
#$path_scripts/Function_calls/subsampling.sh $data/PBMCMosaic_Reference $output 300 \
#102 $path_scripts

#qsub -N abundance -l m_mem_free=32G -l h_rt=03:00:00 \
#$path_scripts/Function_calls/subsampling.sh $data/PBMC10x_Reference $output 20 \ "1022,1373,2418,3090" $path_scripts 
    
#qsub -N abundance -l m_mem_free=32G -l h_rt=03:00:00 \
#$path_scripts/Function_calls/subsampling.sh $data/PBMCMosaic_Reference $output 20 \ "2274,3763,551,4321,6611,209,890,818,102 " $path_scripts 


qsub -N runtime -l m_mem_free=16G -l h_rt=01:00:00 \
$path_scripts/Function_calls/subsampling.sh $data/PBMCMosaic_Reference $output 10 \ "38,100,500,1000,2000,3000" $path_scripts 


