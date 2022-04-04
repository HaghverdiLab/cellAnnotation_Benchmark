
#!/bin/sh


source /home/cmoelbe/bin/anaconda3/bin/activate ItClust


start=`date +%s`
python3 $5/ItClust/reproduce_itClust.py $1/$4 $2/$4 $3/$4
end=`date +%s`

echo ItClust $4 `expr $end - $start` >> /fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection/Results/runtime.txt

