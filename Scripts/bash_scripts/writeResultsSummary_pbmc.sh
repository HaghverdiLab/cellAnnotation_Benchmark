#!/bin/sh

path=/fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection

path_scripts="$path/Scripts/Analysis/"
folder="$path/Data/results"


celltypeFile="$path/Results/celltypeAccuracies_PBMC.csv"
umapFile="$path/Results/umap_PBMC.csv"

source /home/cmoelbe/bin/anaconda3/bin/activate benchmarkpaper 


echo "---------------------- Seurat ----------------------"
#for tag in $folder/Seurat/Mo* ; do
    #python3 $path_scripts/ResultsToTable.py $tag  $celltypeFile Seurat
#    Rscript $path_scripts/getUMAPfile.R Seurat $tag $umapFile 
#done


echo "---------------------- SCN ----------------------"
for tag in $folder/SCN/Mo* ; do
    python3 $path_scripts/ResultsToTable.py $tag  $celltypeFile SCN
    Rscript $path_scripts/getUMAPfile.R SCN $tag $umapFile 
done

echo "---------------------- CellID ----------------------"
#for tag in $folder/CellID/Mo* ; do
   #python3 $path_scripts/ResultsToTable.py $tag  $celltypeFile CellID
#   Rscript $path_scripts/getUMAPfile.R CellID $tag $umapFile CellID
#done

echo "---------------------- MLP ----------------------"
#for tag in $folder/MLP/Mo*; do
   #python3 $path_scripts/ResultsToTable.py  $tag  $celltypeFile MLP
#   Rscript $path_scripts/getUMAPfile.R MLP $tag $umapFile  $path/Data/processed/PBMC_mono/meta_test.txt 
#done


echo "---------------------- ItClust ----------------------"
#for tag in $folder/ItClust/Mosaic* ; do
  #   Rscript $path_scripts/itClust_transformResults.R $tag
  #   Rscript $path_scripts/getUMAPfile.R ItClust $tag $umapFile
 #    echo $tag
#     python3 $path_scripts/ResultsToTable.py $tag $celltypeFile ItClust
#done


#echo "---------------------- SingleR ----------------------"
#for tag in $folder/singleR/Mo* ; do
#    python3 $path_scripts/ResultsToTable.py $tag  $celltypeFile SingleR
#    Rscript $path_scripts/getUMAPfile.R SingleR $tag $umapFile 
#done
