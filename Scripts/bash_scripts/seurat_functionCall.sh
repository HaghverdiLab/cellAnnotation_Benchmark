#!/bin/sh
source  /home/cmoelbe/bin/anaconda3/bin/activate benchmarkpaper
start=`date +%s`

Rscript $5/seurat.R $1/$2 $3 $4/$2

end=`date +%s`

#echo Seurat $2 `expr $end - $start` >> /fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection/Results/runtime.txt



