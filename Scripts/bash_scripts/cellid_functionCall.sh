#!/bin/sh
source  /home/cmoelbe/bin/anaconda3/bin/activate benchmarkpaper
start=`date +%s`

Rscript $5/cellid.R $1 $2 $3 $4

end=`date +%s`

#echo CellID $2 `expr $end - $start` >> /fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection/Results/runtime.txt