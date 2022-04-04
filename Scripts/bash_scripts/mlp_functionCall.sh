#!/bin/sh

source  /home/cmoelbe/bin/anaconda3/bin/activate scPotter

start=`date +%s`
time python3 $4/run_save.py \
                -i $2/ -tag $1 \
                --oo $3/$1/   \
                --classifier MLP \
                --n_hidden_FC 32 \
                --n_hidden_FC2 16 \
                --dropout 0.1 \
                --infer_graph False \
                --seed 0 \
                --K 2 \
                --epochs 30
end=`date +%s`

echo MLP $1 `expr $end - $start` >> /fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection/Results/runtime.txt





