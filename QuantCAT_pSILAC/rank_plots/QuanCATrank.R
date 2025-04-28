library(ggplot2)
library(tidyverse)
library(ggrepel)

neuriteRNAsredo <- read.csv("/Users/ryanvh/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/TTU-niacard02201645/RBPkd_boyden_manuscript/FromDTI/Ziyi/QuantCAT/pSILAC/WT-Soma-Neurite_neuriteenriched.csv")
neuriteQH_rank <- read.csv("/Users/ryanvh/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/TTU-niacard02201645/RBPkd_boyden_manuscript/FromDTI/Ziyi/QuantCAT/QuantCAT/QH_rank.csv")
neuriteQHf_rank <- read.csv("/Users/ryanvh/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/TTU-niacard02201645/RBPkd_boyden_manuscript/FromDTI/Ziyi/QuantCAT/QuantCAT/QHf_rank.csv")
neuriteQHi_rank <- read.csv("/Users/ryanvh/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/TTU-niacard02201645/RBPkd_boyden_manuscript/FromDTI/Ziyi/QuantCAT/QuantCAT/QHi_rank.csv")
selectednamesQ <- read.csv("/Users/ryanvh/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/TTU-niacard02201645/RBPkd_boyden_manuscript/FromDTI/Ziyi/QuantCAT/QuantCAT/selectedgenes.csv")

ggplot(data = neuriteQH_rank, aes(x=Rank, y=log2(row_mean.QH_H))) +
  geom_point(color="#969696") +
  geom_point(data=filter(neuriteQH_rank, Genes %in% neuriteRNAsredo$'SYMBOL'), color = "#7F3F98") +
  geom_label_repel(
    data=neuriteQH_rank %>% 
      filter(Genes %in% neuriteRNAsredo$`SYMBOL`), 
    aes(label=Genes), color="#7F3F98", max.overlaps = Inf) +
  theme_classic()

ggplot(data = neuriteQH_rank, aes(x=Rank, y=log2(row_mean.QH_H))) +
  geom_point(color="#969696") +
  geom_point(data=filter(neuriteQH_rank, Genes %in% neuriteRNAsredo$'SYMBOL'), color = "#7F3F98") +
  geom_label_repel(
    data=neuriteQH_rank %>% 
      filter(Genes %in% selectednamesQ$`SYMBOL`), 
    aes(label=Genes), color="#7F3F98", max.overlaps = Inf) +
  theme_classic()
ggsave('QuanCAT_QH.pdf',width = 6,height = 5)

ggplot(data = neuriteQHf_rank, aes(x=Rank, y=log2(row_mean.QHf_H))) +
  geom_point(color="#969696") +
  geom_point(data=filter(neuriteQHf_rank, Genes %in% neuriteRNAsredo$'SYMBOL'), color = "#7F3F98") +
  geom_label_repel(
    data=neuriteQHf_rank %>% 
      filter(Genes %in% neuriteRNAsredo$`SYMBOL`), 
    aes(label=Genes), color="#7F3F98", max.overlaps = Inf) +
  theme_classic()
ggsave('QuanCAT_QHf.pdf',width = 6,height = 5)

ggplot(data = neuriteQHi_rank, aes(x=Rank, y=log2(row_mean.QHi_H))) +
  geom_point(color="#969696") +
  geom_point(data=filter(neuriteQHi_rank, Genes %in% neuriteRNAsredo$'SYMBOL'), color = "#7F3F98") +
  geom_label_repel(
    data=neuriteQHi_rank %>% 
      filter(Genes %in% neuriteRNAsredo$`SYMBOL`), 
    aes(label=Genes), color="#7F3F98", max.overlaps = Inf) +
  theme_classic()
ggsave('QuanCAT_QHi.pdf',width = 6,height = 5)
