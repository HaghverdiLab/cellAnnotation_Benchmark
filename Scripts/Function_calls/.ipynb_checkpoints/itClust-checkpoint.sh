#!/bin/sh



Filelist=($2$(ls -d $2/$1* | xargs -n1  basename))

source  /home/cmoelbe/bin/anaconda3/bin/activate r-environment
Rscript Scripts/Methods/prep_ItClust.R  ${Filelist[${SGE_TASK_ID}-1]} $2 $3 $4

source /home/cmoelbe/bin/anaconda3/bin/activate ItClust

start=`date +%s`
python3 Scripts/Methods/itClust.py $4/${Filelist[${SGE_TASK_ID}-1]} $4/${Filelist[${SGE_TASK_ID}-1]} $4/${Filelist[${SGE_TASK_ID}-1]}
end=`date +%s`

echo ItClust ${Filelist[${SGE_TASK_ID}-1]} `expr $end - $start` >> /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Results/Files/runtime.txt

Rscript Scripts/itClust_transformResults.R $4/${Filelist[${SGE_TASK_ID}-1]}