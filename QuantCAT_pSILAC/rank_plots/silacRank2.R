library(ggplot2)
library(tidyverse)
library(ggrepel)

pSILAC <- read.csv("/Users/ryanvh/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/TTU-niacard02201645/RBPkd_boyden_manuscript/FromDTI/Ziyi/QuantCAT/pSILAC/mean_value.csv")
neuritepSILAC <- read.csv("/Users/ryanvh/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/TTU-niacard02201645/RBPkd_boyden_manuscript/FromDTI/Ziyi/QuantCAT/pSILAC/mean_value_neuritedetected.csv")
neuriteRNAs <- read.csv("/Users/ryanvh/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/TTU-niacard02201645/RBPkd_boyden_manuscript/FromDTI/Ziyi/QuantCAT/pSILAC/lfc05_zone_WT.csv")
neuriteRNAsredo <- read.csv("/Users/ryanvh/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/TTU-niacard02201645/RBPkd_boyden_manuscript/FromDTI/Ziyi/QuantCAT/pSILAC/WT-Soma-Neurite_neuriteenriched.csv")
selectednames <- read.csv("/Users/ryanvh/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/TTU-niacard02201645/RBPkd_boyden_manuscript/FromDTI/Ziyi/QuantCAT/pSILAC/selectednames.csv")

ggplot(data = neuritepSILAC, aes(x=rank, y=log2(row_mean.HN_H))) +
  geom_point(color="#969696") +
  geom_point(data=filter(neuritepSILAC, PG.Genes %in% neuriteRNAsredo$'SYMBOL'), color = "#7F3F98") +
  geom_label_repel(
    data=neuritepSILAC %>% 
      filter(PG.Genes %in% neuriteRNAsredo$`SYMBOL`) %>%
      filter(rank < 80), 
    aes(label=PG.Genes), color="#7F3F98", max.overlaps = Inf) +
  geom_label_repel(
    data=neuritepSILAC %>% 
      filter(PG.Genes %in% neuriteRNAsredo$`SYMBOL`) %>%
      filter(rank > 2033), 
    aes(label=PG.Genes), color="#7F3F98", max.overlaps = Inf) +
  theme_classic()
ggsave('pSILACrankredo.pdf',width = 6,height = 5)


ggplot(data = neuritepSILAC, aes(x=rank, y=log2(row_mean.HN_H))) +
  geom_point(color="#969696") +
  geom_point(data=filter(neuritepSILAC, PG.Genes %in% neuriteRNAsredo$'SYMBOL'), color = "#7F3F98") +
  geom_label_repel(
    data=neuritepSILAC %>% 
      filter(PG.Genes %in% selectednames$`SYMBOL`),
    aes(label=PG.Genes), color="#7F3F98", max.overlaps = Inf) +
  theme_classic()
ggsave('pSILACrankredo.pdf',width = 6,height = 5)
