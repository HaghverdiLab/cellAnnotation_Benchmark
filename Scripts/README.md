# Scripts used in the project
The scripts are used in individual conda enviroments based on the installation information found on their representing websites. The GRN construction scripts use the same environment as Seurat. Since they use a big overlap in functions.

## Analysis
Scripts used to analyse and visualize the predictions made with the compared methods.

## Data preparation
Scripts used to prepare the PBMC and the Cross Species data set for the analysis.

## ScPotter
Scripts to use scPotter. This analysis focusses on understanding the benefits and weaknesses of scPotter and to find possible improvements. More information about the underlying method can be found here https://github.com/EliHei2/scPotter. In order to run scPotter, scPotte must be set up on the machine. Here, we did this in side a conda environment which is called in the functionCall script.


## ItClust
Scripts to use the ItClust. More information about ItClust can be found here https://github.com/jianhuupenn/ItClust

## Seurat
Scripts to use Seurat. More information about Seurat can be found here https://satijalab.org/seurat/

## MLP
MLP is a method implemented in the scPotter package as a controll method not relying on the gene regulatory network
