#!/bin/sh
#$ -N ItClust_$1
#$ -l m_mem_free=12G

source  /home/cmoelbe/bin/anaconda3/bin/activate scPotterPaper

Rscript $5/ItClust/getItClust.R $1 $2 $3 $4
