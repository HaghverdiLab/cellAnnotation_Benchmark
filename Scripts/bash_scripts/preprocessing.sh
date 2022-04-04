#!/bin/sh
path=/fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection
path_scripts="$path/Scripts"
path_data="$path/Data/subsets/"
output="$path/Data/results/"
echo "Starting"
test="$path/Data/processed/PBMC_mono"

for x in $(ls -d $path_data/Mo*| xargs -n1  basename); do    
         echo $x
         #qsub -N $x -l m_mem_free=32G -l h_rt=0:30:00 -o itclust_out_$x  -e itclust_error_$x \
         #$path_scripts/bash_scripts/itClust_getFiles.sh $x $path_data $test $output/ItClust/ $path_scripts/methods
         
         #qsub -N $x -l m_mem_free=32G -l h_rt=0:30:00  -o mlp_out_$x  -e mlp_error_$x \
         #$path_scripts/bash_scripts/mlp_prepare_functionCall.sh $path_data $x $test $output/MLP $path_scripts/methods
         
         qsub -N $x -l m_mem_free=32G -l h_rt=0:30:00  -pe smp 4 -o seurat_out_$x  -e seurat_error_$x \
         $path_scripts/bash_scripts/seurat_functionCall.sh $path_data $x $test $output/Seurat $path_scripts/methods
         
         #qsub -N $x -l m_mem_free=32G -l h_rt=0:30:00  -pe smp 4 -o cellid_out_$x  -e cellid_error_$x \
         #$path_scripts/bash_scripts/cellid_functionCall.sh $path_data $x $test $output/CellID $path_scripts/methods
         
         #qsub -N $x -l m_mem_free=32G -l h_rt=0:10:00  -pe smp 4 -o scn_out_$x  -e scn_error_$x \
         #$path_scripts/bash_scripts/scn_functionCall.sh $path_data $x $test $output/SCN $path_scripts/methods
         
         #qsub -N $x -l m_mem_free=32G -l h_rt=0:30:00  -pe smp 4 -o singler_out_$x  -e singler_error_$x \
          #$path_scripts/bash_scripts/singleR_functionCall.sh $path_data $x $test $output/singleR $path_scripts/methods
done 

test="$path/Data/processed/CrossSpecies/"
for x in $(ls -d $path_data/CrossSpecies_228* | xargs -n1  basename); do    
        echo $x
         #qsub -N $x -l m_mem_free=32G -l h_rt=0:30:00  -o itclust_out_$x  -e itclust_error_$x \
         #$path_scripts/bash_scripts/itClust_getFiles.sh $x $path_data $test $output/ItClust/ $path_scripts/methods
         
         #qsub -N $x -l m_mem_free=32G -l h_rt=0:30:00  -o mlp_out_$x  -e mlp_error_$x \
         #$path_scripts/bash_scripts/mlp_prepare_functionCall.sh $path_data $x $test $output/MLP $path_scripts/methods
         
         qsub -N $x -l m_mem_free=32G -l h_rt=0:30:00  -pe smp 4 -o seurat_out_$x  -e seurat_error_$x \
         $path_scripts/bash_scripts/seurat_functionCall.sh $path_data/ $x $test  $output/Seurat $path_scripts/methods
         
         #qsub -N $x -l m_mem_free=32G -l h_rt=0:30:00  -pe smp 4 -o cellid_out_$x  -e cellid_error_$x \
         #$path_scripts/bash_scripts/cellid_functionCall.sh $path_data $x $test $output/CellID $path_scripts/methods
         
         #qsub -N $x -l m_mem_free=32G -l h_rt=0:10:00  -pe smp 4 -o scn_out_$x  -e scn_error_$x \
         #$path_scripts/bash_scripts/scn_functionCall.sh $path_data $x $test $output/SCN $path_scripts/methods
         
         #qsub -N $x -l m_mem_free=32G -l h_rt=0:30:00  -pe smp 4 -o singler_out_$x  -e singler_error_$x \
         #$path_scripts/bash_scripts/singleR_functionCall.sh $path_data $x $test $output/singleR $path_scripts/methods
done 


