devtools::document('../../../lib/lcdbwf')
devtools::load_all('../../../lib/lcdbwf')
library(dplyr)
library(AnnotationHub)
library(BiocParallel)
library(cowplot)
library(DESeq2)
library(dplyr)
library(DT)
library(genefilter)
library(ggplot2)
library(gridExtra)
library(plotly)
library(purrr)
library(readr)
library(reshape)
library(tibble)
library(tximport)
library(tidyr)
# NOTE: requires transite env because of different version of DESeq for plotting

## Draws the PCA
## first_PC numeric value of the x-axis PC to draw
## second_PC numeric value of the y-axis PC to draw
## group_name: A name for the settings of the plot. Appears in the title and filename
## my_vsd: a subset vsd for plotting
make_PCA <- function(first_PC, second_PC, group_name, my_vsd){

  # PCA without correction
  pca_matrix<- as.data.frame(assay(my_vsd)) %>%
      base::as.matrix() %>%
      t()
  sample_pca <- prcomp(pca_matrix)

  pc_scores <- sample_pca$x %>% as_tibble(rownames="sample")

  # Change variables in order to generate PCA plots
  PC1 <- first_PC
  PC2 <- second_PC
  PC_first <- paste0('PC', PC1)
  PC_second <- paste0('PC', PC2)

  mat <- DESeq2::plotPCA(my_vsd, c("group", 'zone', 'fillcolor', 'targeting'), pcsToUse=PC1:PC2, returnData=TRUE)
  pv <- attr(mat, 'percentVar')
  
  p <- ggplot(data=mat,aes_string(x = PC_first, y = PC_second, color='group.1', shape='zone', fill='fillcolor')) + 
       geom_point(size=3) +
       scale_shape_manual(values=c('N'=21, 'S'=22)) +
       scale_colour_manual(values=c('FUSkd'='#FF89D7','TDP43kd' ='#D682FF', 'WT'='black', 'HNRNPA1kd'='#0095FF')) +
       scale_fill_manual(values=c('#FF89D7'='#FF89D7','#D682FF' ='#D682FF', 'black'='black', '#0095FF'='#0095FF'), na.value='transparent') +
       theme_minimal() + xlab(paste0(PC_first,': ', round(pv[1]*100), '% variance')) +
       guides(fill='none',
             shape = guide_legend(override.aes = list(fill = "black"))) + 
       ylab(paste0(PC_second,': ', round(pv[2]*100), '% variance')) + coord_fixed() +
       ggtitle(paste0("PCA: ", group_name," empty shapes are non-targeting samples", " , ", PC_first, " and ", PC_second))

  ggsave(paste0("PCA_plots/cory_style_final/", group_name, "_", "zone_trgt_", "pca_",PC_first, "_", PC_second,".pdf"), plot=p)
}

config <- lcdbwf:::load_config('config.yaml')

parallel <- config$parallel$parallel
if (config$parallel$parallel){
    register(MulticoreParam(config$parallel$cores))
}


colData <- read.table(config$main$sampletable, sep='\t', header=TRUE, stringsAsFactors=FALSE)

# lcdb-wf requires that the first column of the sampletable contain sample
# names. Use that to set the rownames.
rownames(colData) <- colData[,1]

colData$targeting <- case_match(colData$guide_group, 
                                     c('sg100', 'sg1126') ~ 'nontargeting',
                                     c('sg200', 'sg1152', 'sg1128') ~ 'targeting',
                                     c('WT') ~ 'WT',
                                     )
colData$targeting <- factor(colData$targeting, levels=c('nontargeting', 'targeting', 'WT'))

# The soma is the reference lvl since it includes everything
colData$zone <- factor(colData$zone, levels=c('S', 'N'))
colData$pair <- factor(colData$pair)

colData_init <- colData

dds_initial <- lcdbwf:::make_dds(
  list(sampletable=colData_init, design=~pair + zone, subset_counts=FALSE),
  config=config,
  parallel=config$parallel$parallel
)
vsd <- varianceStabilizingTransformation(dds_initial, blind=FALSE)

vsd$fillcolor[(vsd$group == 'FUSkd')] <- '#FF89D7'
vsd$fillcolor[(vsd$group == 'HNRNPA1kd')] <- '#0095FF'
vsd$fillcolor[(vsd$group == 'TDP43kd')] <- '#D682FF'
vsd$fillcolor[(vsd$group == 'WT')] <- 'black'
vsd$fillcolor[(vsd$targeting == 'nontargeting')] <- NA

vsd$fillcolor <- factor(vsd$fillcolor)
vsd$group <- factor(vsd$group)


# FUS set
current_vsd <- vsd[,grepl('FUSkd',vsd$group)]

make_PCA(1, 2, "FUSkd", current_vsd)
make_PCA(2, 3, "FUSkd", current_vsd)
make_PCA(3, 4, "FUSkd", current_vsd)

# HNRNPA1kd set
current_vsd <- vsd[,grepl('HNRNPA1kd',vsd$group)]

make_PCA(1, 2, "HNRNPA1kd", current_vsd)
make_PCA(2, 3, "HNRNPA1kd", current_vsd)
make_PCA(3, 4, "HNRNPA1kd", current_vsd)

# HNRNPA1kd set
current_vsd <- vsd[,grepl('TDP43kd',vsd$group)]
make_PCA(1, 2, "TDP43kd", current_vsd)
make_PCA(2, 3, "TDP43kd", current_vsd)
make_PCA(3, 4, "TDP43kd", current_vsd)

# No-WT
current_vsd <- vsd[,vsd$group !='WT']

make_PCA(1, 2, "all_except_WT", current_vsd)
make_PCA(2, 3, "all_except_WT", current_vsd)
make_PCA(3, 4, "all_except_WT", current_vsd)

# All samples
current_vsd <- vsd
make_PCA(1, 2, "all_samples", current_vsd)
make_PCA(2, 3, "all_samples", current_vsd)
make_PCA(3, 4, "all_samples", current_vsd)
