folder=/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/
folder_data=$folder"/Data/"
meta=$folder_data"/Fulldata/PBMC_Query/meta.csv"
name=PBMC10x
folder_results=$folder"/Results/Files"

#Rscript $path/Notebooks/get_result_file.r $folder_data"/Predictions/" $meta $name $folder_results"/result_general.csv" #$folder_results"/results_long.csv"


#Rscript $path/Notebooks/get_confidence_scores.r $folder_data"/Predictions/" $name 38,100,100

Rscript $path/Notebooks/get_umap_file.r $name"_100_" PBMCMosaicBalanced_102_ $folder_data"/Predictions/" $meta /fast$folder_results"/result_general.csv" $folder_results"/results_curated_PBMC.csv"