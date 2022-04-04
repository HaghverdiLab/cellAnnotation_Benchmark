#!/bin/sh
source /home/cmoelbe/bin/anaconda3/bin/activate scPotterPaper

Rscript /fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection/Scripts/Data_Preparation/create_PBMC_sets.R 
Rscript /fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection/Scripts/Data_Preparation/create_crossSpecies_sets.r 
