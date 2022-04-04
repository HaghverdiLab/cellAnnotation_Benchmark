#!/bin/sh
path=/fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection
path_scripts="$path/Scripts"
output="$path/Data/results"



#for tag in $(ls -d $output/ItClust/Mosaic_3000* | xargs -n1  basename) ; do
#   qsub -N $tag -l m_mem_free=32G  -pe smp 4 -o /home/cmoelbe/out/itClust_$tag -e /home/cmoelbe/out/itClust_error_$tag\
#   $path_scripts/bash_scripts/itClust_functionCall.sh $output/ItClust/  $output/ItClust/ \
#   $output/ItClust/  $tag/ $path_scripts
#done

for tag in $(ls -d $output/MLP/* | xargs -n1  basename) ; do
    qsub -N $tag -l m_mem_free=32G  -pe smp 4 -o $output/MLP/$tag/out_$tag -e $output/MLP/$tag/error_$tag \
    $path_scripts/bash_scripts/mlp_functionCall.sh $tag  $output/MLP/ $output/MLP $path_scripts/scPotter/
done
