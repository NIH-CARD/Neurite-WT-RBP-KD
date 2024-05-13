library(dplyr)
library(transite)

names <- c("sg200_TDP43kd", "sg100_TDP43kd", "sg1126_HNRNPA1kd", 
        "sg100_FUSkd", "sg1152_FUSkd", "sg1128_HNRNPA1kd")

base_folder <- 'transite_logfc_results'

if (!file.exists(base_folder)){dir.create(base_folder)}
for (name in names){
    object_path <- paste0(name,".rds")
    
    result <- readRDS(object_path)
    write.table(result$spectrum_info_df, paste0(base_folder, '/', name, "_result.csv"), quote=FALSE, sep=';')
    
    # create plots
    i <- 0
    
    for (plotx in result$spectrum_plots){
        i <- i+1
        if (!(file.exists(file.path(base_folder, paste0(name,"_plots"))))){
                  dir.create(file.path(base_folder, paste0(name, "_plots")))
                }
        pdf(paste0(base_folder, '/', name, "_plots/", i, '.pdf'))
        plot(plotx)
        dev.off()
    }
}
