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
# requires transite env because of different version of DESeq for plotting

obj <- readRDS('combined.Rds')
res_list <- obj$res_list
dds_list <- obj$dds_list

parallel <- config$parallel$parallel
if (config$parallel$parallel){
    register(MulticoreParam(config$parallel$cores))
}

config <- lcdbwf:::load_config('config.yaml')

colData <- read.table(config$main$sampletable, sep='\t', header=TRUE, stringsAsFactors=FALSE)

# lcdb-wf requires that the first column of the sampletable contains sample
# names. Use that to set the rownames.
rownames(colData) <- colData[,1]

colData$targeting <- case_match(colData$guide_group, 
                                     c('sg100', 'sg1126') ~ 'nontargeting',
                                     c('sg200', 'sg1152', 'sg1128') ~ 'targeting',
                                     c('WT') ~ 'WT',
                                     )
colData$targeting <- factor(colData$targeting, levels=c('nontargeting', 'targeting', 'WT'))
colData$pair <- factor(colData$pair)
# The soma is the reference lvl since it includes everything
colData$zone <- factor(colData$zone, levels=c('S', 'N'))
colData <- colData %>% mutate(batch = factor(case_when(group == 'TDP43kd' ~ 'batch1',
                                                             group == 'FUSkd' ~ 'batch2',
                                                             group == 'HNRNPA1kd' ~ 'batch3',
                                                             group == 'WT' ~ 'WT')))

# check if 'rld_list' exists in the object, if not sets to NULL
if(!'rld_list' %in% names(obj)){
  rld_list <- NULL
} else {
  rld_list <- obj$rld_list
}

colData_init <- colData[colData$'group'=='FUSkd',]
colData_init <- colData

dds_initial <- lcdbwf:::make_dds(
  list(sampletable=colData_init, design=~1, subset_counts=TRUE),
  config=config,
  parallel=config$parallel$parallel
)
vsd <- varianceStabilizingTransformation(dds_initial, blind=TRUE)


make_PCA <- function(first_PC, second_PC, gene, var,terms, color_vector){


  vsd_no_wt <- vsd[,grepl(terms,vsd$group)]


  #PCA without correction
  pca_matrix<- as.data.frame(assay(vsd_no_wt)) %>%
      base::as.matrix() %>%
      t()
  sample_pca <- prcomp(pca_matrix)

  # HAven't satandarzied yet
  pc_scores <- sample_pca$x %>% as_tibble(rownames="sample")

  # Change variables in order to generate PCA plots
  variable <- var
  PC1 <- first_PC
  PC2 <- second_PC
  first_PC <- paste0('PC', PC1)
  second_PC <- paste0('PC', PC2)

# Set shapes
  vsd_no_wt$shapes[(vsd_no_wt$targeting %in% c('targeting', 'WT')) & vsd_no_wt$zone == 'N'] <- '1'
  vsd_no_wt$shapes[(vsd_no_wt$targeting == 'targeting') & vsd_no_wt$zone == 'S'] <- '2'
  vsd_no_wt$shapes[(vsd_no_wt$targeting == 'nontargeting') & vsd_no_wt$zone == 'N'] <- '3'
  vsd_no_wt$shapes[(vsd_no_wt$targeting %in% c('nontargeting','WT')) & vsd_no_wt$zone == 'S'] <- '4'

  zone_var <- if(variable == 'group') 'group.1' else variable # group is similar to the internal var

  pc_metadata <- inner_join(pc_scores, as.tibble(colData(vsd_no_wt)), join_by(sample == sample))
  mat <- DESeq2::plotPCA(vsd_no_wt, c(variable, 'targeting'), pcsToUse=PC1:PC2, returnData=TRUE)
  pv <- attr(mat, 'percentVar')
  
  p <- ggplot(data=mat,aes_string(x = first_PC, y = second_PC, color="targeting", shape=zone_var)) +
      scale_shape_manual(values=c('1'=1, '2'=16, '3'=0,'4'=15)) +
      scale_colour_manual(values=color_vector) +
      theme_bw()+ 
  geom_point(size=3) + xlab(paste0(first_PC,': ', round(pv[1]*100), '% variance')) +
            ylab(paste0(second_PC,': ', round(pv[2]*100), '% variance')) + coord_fixed() +
            ggtitle(paste0(variable, " , ", first_PC, " and ", second_PC, " ", gene))
  ggsave(paste0("PCA_plots/final_plots/", gene, "_", variable, "_pca_",first_PC, "_", second_PC,".pdf"), plot=p)
}

make_PCA(1, 2, "FUSkd", "shapes", "FUS", color_vector=c("targeting"="#FF89D7", "nontargeting"="black"))
make_PCA(1, 2, "HNRNPA1kd", "shapes", "HNRN", color_vector=c("targeting"="#0095FF", "nontargeting"="black"))
make_PCA(1, 2, "TDP43kd", "shapes", "TDP", color_vector=c("targeting"="#D682FF", "nontargeting"="black"))
make_PCA(1, 2, "WT", "shapes", "WT", color_vector=c("WT"="black"))

make_PCA(2, 3, "FUSkd", "shapes", "FUS", color_vector=c("targeting"="#FF89D7", "nontargeting"="black"))
make_PCA(2, 3, "HNRNPA1kd", "shapes", "HNRN", color_vector=c("targeting"="#0095FF", "nontargeting"="black"))
make_PCA(2, 3, "TDP43kd", "shapes", "TDP", color_vector=c("targeting"="#D682FF", "nontargeting"="black"))
make_PCA(2, 3, "WT", "shapes", "WT", color_vector=c("WT"="black"))


make_PCA(3, 4, "FUSkd", "shapes", "FUS", color_vector=c("targeting"="#FF89D7", "nontargeting"="black"))
make_PCA(3, 4, "HNRNPA1kd", "shapes", "HNRN", color_vector=c("targeting"="#0095FF", "nontargeting"="black"))
make_PCA(3, 4, "TDP43kd", "shapes", "TDP", color_vector=c("targeting"="#D682FF", "nontargeting"="black"))
make_PCA(3, 4, "WT", "shapes", "WT", color_vector=c("WT"="black"))

  vsd_no_wt <- vsd[,vsd$group !='WT']


  #PCA without correction
  pca_matrix<- as.data.frame(assay(vsd_no_wt)) %>%
      base::as.matrix() %>%
      t()
  sample_pca <- prcomp(pca_matrix)

  # HAven't satandarzied yet
  pc_scores <- sample_pca$x %>% as_tibble(rownames="sample")

    vsd_no_wt$shapes[(vsd_no_wt$targeting %in% c('targeting', 'WT')) & vsd_no_wt$zone == 'N'] <- '1'
  vsd_no_wt$shapes[(vsd_no_wt$targeting == 'targeting') & vsd_no_wt$zone == 'S'] <- '2'
  vsd_no_wt$shapes[(vsd_no_wt$targeting == 'nontargeting') & vsd_no_wt$zone == 'N'] <- '3'
  vsd_no_wt$shapes[(vsd_no_wt$targeting %in% c('nontargeting','WT')) & vsd_no_wt$zone == 'S'] <- '4'

  # Change variables in order to generate PCA plots
  variable <- 'zone'
  PC1 <- 2
  PC2 <- 3
  first_PC <- paste0('PC', PC1)
  second_PC <- paste0('PC', PC2)

  vsd_no_wt$color[(vsd_no_wt$targeting == 'nontargeting')] <- 'black'
  vsd_no_wt$color[(vsd_no_wt$targeting == 'targeting') & vsd_no_wt$group == 'TDP43kd'] <- '#D682FF'
  vsd_no_wt$color[(vsd_no_wt$targeting == 'targeting') & vsd_no_wt$group == 'FUSkd'] <- '#FF89D7'
  vsd_no_wt$color[(vsd_no_wt$targeting == 'targeting') & vsd_no_wt$group == 'HNRNPA1kd'] <- '#0095FF'

  vsd_no_wt$shapes[(vsd_no_wt$targeting %in% c('targeting', 'WT')) & vsd_no_wt$zone == 'N'] <- '1'
  vsd_no_wt$shapes[(vsd_no_wt$targeting == 'targeting') & vsd_no_wt$zone == 'S'] <- '2'
  vsd_no_wt$shapes[(vsd_no_wt$targeting == 'nontargeting') & vsd_no_wt$zone == 'N'] <- '3'
  vsd_no_wt$shapes[(vsd_no_wt$targeting %in% c('nontargeting','WT')) & vsd_no_wt$zone == 'S'] <- '4'

  pc_metadata <- inner_join(pc_scores, as.tibble(colData(vsd_no_wt)), join_by(sample == sample))
  mat <- DESeq2::plotPCA(vsd_no_wt, c('color', 'shapes'), pcsToUse=PC1:PC2, returnData=TRUE)
  pv <- attr(mat, 'percentVar')

    p <- ggplot(data=mat,aes_string(x = first_PC, y = second_PC, color = 'color', shape='shapes')) +
        scale_shape_manual(values=c('1'=1, '2'=16, '3'=0,'4'=15)) +
      scale_colour_manual(values=unique(vsd_no_wt$color)) +
      theme_bw()+ 
    geom_point(size=3) + xlab(paste0(first_PC,': ', round(pv[1]*100), '% variance')) +
              ylab(paste0(second_PC,': ', round(pv[2]*100), '% variance')) + coord_fixed() +
              ggtitle(paste0(variable, " , ", first_PC, " and ", second_PC, " all non-WT samples")) + 
              theme(panel.background = element_blank())
    ggsave(paste0("PCA_plots/final_plots/", "pca_",first_PC, "_", second_PC,"_all_no_WT.pdf"), plot=p)



##### With WT
  vsd_no_wt <- vsd


  #PCA without correction
  pca_matrix<- as.data.frame(assay(vsd_no_wt)) %>%
      base::as.matrix() %>%
      t()
  sample_pca <- prcomp(pca_matrix)

  # HAven't satandarzied yet
  pc_scores <- sample_pca$x %>% as_tibble(rownames="sample")

    vsd_no_wt$shapes[(vsd_no_wt$targeting %in% c('targeting', 'WT')) & vsd_no_wt$zone == 'N'] <- '1'
  vsd_no_wt$shapes[(vsd_no_wt$targeting == 'targeting') & vsd_no_wt$zone == 'S'] <- '2'
  vsd_no_wt$shapes[(vsd_no_wt$targeting == 'nontargeting') & vsd_no_wt$zone == 'N'] <- '3'
  vsd_no_wt$shapes[(vsd_no_wt$targeting %in% c('nontargeting','WT')) & vsd_no_wt$zone == 'S'] <- '4'

  # Change variables in order to generate PCA plots
  variable <- 'zone'
  PC1 <- 2
  PC2 <- 3
  first_PC <- paste0('PC', PC1)
  second_PC <- paste0('PC', PC2)

  vsd_no_wt$color[(vsd_no_wt$targeting == 'nontargeting')] <- 'black'
  vsd_no_wt$color[(vsd_no_wt$targeting == 'targeting') & vsd_no_wt$group == 'TDP43kd'] <- '#D682FF'
  vsd_no_wt$color[(vsd_no_wt$targeting == 'targeting') & vsd_no_wt$group == 'FUSkd'] <- '#FF89D7'
  vsd_no_wt$color[(vsd_no_wt$targeting == 'targeting') & vsd_no_wt$group == 'HNRNPA1kd'] <- '#0095FF'

  vsd_no_wt$shapes[(vsd_no_wt$targeting %in% c('targeting', 'WT')) & vsd_no_wt$zone == 'N'] <- '1'
  vsd_no_wt$shapes[(vsd_no_wt$targeting == 'targeting') & vsd_no_wt$zone == 'S'] <- '2'
  vsd_no_wt$shapes[(vsd_no_wt$targeting == 'nontargeting') & vsd_no_wt$zone == 'N'] <- '3'
  vsd_no_wt$shapes[(vsd_no_wt$targeting %in% c('nontargeting','WT')) & vsd_no_wt$zone == 'S'] <- '4'

  pc_metadata <- inner_join(pc_scores, as.tibble(colData(vsd_no_wt)), join_by(sample == sample))
  mat <- DESeq2::plotPCA(vsd_no_wt, c('color', 'shapes'), pcsToUse=PC1:PC2, returnData=TRUE)
  pv <- attr(mat, 'percentVar')

    p <- ggplot(data=mat,aes_string(x = first_PC, y = second_PC, color = 'color', shape='shapes')) +
        scale_shape_manual(values=c('1'=1, '2'=16, '3'=0,'4'=15)) +
      scale_colour_manual(values=unique(vsd_no_wt$color)) +
      theme_bw()+ 
    geom_point(size=3) + xlab(paste0(first_PC,': ', round(pv[1]*100), '% variance')) +
              ylab(paste0(second_PC,': ', round(pv[2]*100), '% variance')) + coord_fixed() +
              ggtitle(paste0(variable, " , ", first_PC, " and ", second_PC, " all")) + 
              theme(panel.background = element_blank())
    ggsave(paste0("PCA_plots/final_plots/", "pca_",first_PC, "_", second_PC,"_all_samples.pdf"), plot=p)