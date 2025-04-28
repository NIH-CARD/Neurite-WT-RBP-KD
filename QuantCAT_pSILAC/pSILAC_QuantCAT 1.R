#########################
#pSALIC and QuantCAT data analysis
#July-2024
#Ziyi Li
#CARD NIH
#########################

setwd('/Users/liz36/Documents/Collaboration/QuantCAT/')
#PACKAGES######################################################################################
package_list = c('ggplot2', 'data.table','reshape2','dplyr',
                 'org.Hs.eg.db','clusterProfiler','limma','ggridges','eulerr',
                 '')
lapply(package_list, require, character.only=TRUE)
source('/Users/liz36/Documents/GitHub/CARD_Protpipe_dev/ProtPipe/src/functions.R')

#DATA INPUT################################
psalic <- fread('SILAC_tSILAC_tSILACReport.csv')
quantcat <- fread('QUANCAT_TSILAC_tSILACReport_update.csv')

#pSALIC#############################
##CLEAN DATA##################################
psalic_h_clean=psalic %>%
  select(1, 2, matches('HN|HS'))%>%
  mutate(across(contains("raw"), ~ ifelse(. < 10, NA, .)))
psalic_h_clean_long=melt_intensity_table(psilac_standardize_format(psalic_h_clean))
plot_silac_pg_intensity(DT.long = psalic_h_clean_long,
                        output_dir = 'pSILAC/QC/HS_HN',
                        height = 4,
                        width = 8
                        )
##Ven plot#################

###data####
HN_L=psalic_h_clean %>%
  select(1, 2, matches('HN.*Channel1'))%>%
  filter(rowSums(is.na(select(., matches("HN.*Channel1")))) < (ncol(select(., matches("HN.*Channel1"))) / 2))%>%
  rowwise() %>%
  mutate(row_mean = mean(c_across(matches("Channel")), na.rm = TRUE))
HS_L=psalic_h_clean %>%
  select(1, 2, matches('HS.*Channel1'))%>%
  filter(rowSums(is.na(select(., matches("HS.*Channel1")))) < (ncol(select(., matches("HS.*Channel1"))) / 2))%>%
  rowwise() %>%
  mutate(row_mean = mean(c_across(matches("Channel")), na.rm = TRUE))
HN_H=psalic_h_clean %>%
  select(1, 2, matches('HN.*Channel3'))%>%
  filter(rowSums(is.na(select(., matches("HN.*Channel3")))) < (ncol(select(., matches("HN.*Channel3"))) / 2))%>%
  rowwise() %>%
  mutate(row_mean = mean(c_across(matches("Channel")), na.rm = TRUE))
HS_H=psalic_h_clean %>%
  select(1, 2, matches('HS.*Channel3'))%>%
  filter(rowSums(is.na(select(., matches("HS.*Channel3")))) < (ncol(select(., matches("HS.*Channel3"))) / 2))%>%
  rowwise() %>%
  mutate(row_mean = mean(c_across(matches("Channel")), na.rm = TRUE))

###plot###########
library(eulerr)
ven_pro=list(
  HN_L=unique(na.omit(HN_L$PG.Genes)),
  HS_L=unique(na.omit(HS_L$PG.Genes)),
  HN_H=unique(na.omit(HN_H$PG.Genes)),
  HS_H=unique(na.omit(HS_H$PG.Genes)))
vd=euler(ven_pro)
pdf("pSILAC/DE/ven_all.pdf",height = 6,width = 7)
par(mar = c(2, 5, 5, 2))
plot(vd, 
     factor_names = TRUE, labels=list(font=2, cex=1.2),
     counts = TRUE,
     key=TRUE,
     cex=1, fills = list(fill = c("#fbb4ae", "#b3cde3", "#ccebc5",'#dfc27d')),
     edges = FALSE,
     quantities = list(type = "counts", "percent"))
dev.off()

ven_pro=list(
  HN_L=unique(na.omit(HN_L$PG.Genes)),
  HN_H=unique(na.omit(HN_H$PG.Genes)))
vd=euler(ven_pro)
pdf("pSILAC/DE/HN.pdf",height = 6,width = 7)
par(mar = c(2, 5, 5, 2))
plot(vd, 
     factor_names = TRUE, labels=list(font=2, cex=1.2),
     counts = TRUE,
     key=TRUE,
     cex=1, fills = list(fill = c("#fbb4ae", "#b3cde3")),
     edges = FALSE,
     quantities = list(type = "counts", "percent"))
dev.off()


ven_pro=list(
  
  HS_L=unique(na.omit(HS_L$PG.Genes)),
 
  HS_H=unique(na.omit(HS_H$PG.Genes)))
vd=euler(ven_pro)
pdf("pSILAC/DE/ven_HS.pdf",height = 6,width = 7)
par(mar = c(2, 5, 5, 2))
plot(vd, 
     factor_names = TRUE, labels=list(font=2, cex=1.2),
     counts = TRUE,
     key=TRUE,
     cex=1, fills = list(fill = c( "#ccebc5",'#dfc27d')),
     edges = FALSE,
     quantities = list(type = "counts", "percent"))
dev.off()


###

###list
# Combine all datasets into one
library(UpSetR)
# Assuming all_data is correctly prepared
# Combine all genes into a single list
all_genes <- unique(c(HN_L$PG.Genes, HS_L$PG.Genes, HN_H$PG.Genes, HS_H$PG.Genes))

# Create a dataframe with genes as rows and datasets as columns
gene_df <- data.frame(
  Gene = all_genes,
  HN_L = as.integer(all_genes %in% HN_L$PG.Genes),
  HS_L = as.integer(all_genes %in% HS_L$PG.Genes),
  HN_H = as.integer(all_genes %in% HN_H$PG.Genes),
  HS_H = as.integer(all_genes %in% HS_H$PG.Genes))

upset(
  all_data,
  sets = c("HN_L", "HS_L", "HN_H", "HS_H"),
  order.by = "freq",
  keep.order = TRUE,
  main.bar.color = "skyblue",
  matrix.color = "lightgray",
  sets.bar.color = "skyblue",
  text.scale = c(1.5, 1, 1.5),
  point.size = 3
)

merged_data <- full_join(HN_L, HS_L, by = "PG.Genes", suffix = c(".HN_L", ".HS_L")) %>%
  full_join(HN_H, by = "PG.Genes") %>%
  full_join(HS_H, by = "PG.Genes", suffix = c(".HN_H", ".HS_H")) %>%
  select(PG.Genes,PG.ProteinGroups, row_mean.HN_L, row_mean.HS_L, row_mean.HN_H, row_mean.HS_H)

write.csv(merged_data,'pSILAC/mean_value.csv')

##DE MA and vocanol #######################
library(limma)
library(ggpubr)

###HN##############
HN=(psalic_h_clean) %>%
  select( matches('HN'))
setrownames(HN, psalic_h_clean$PG.ProteinGroups)
row.names(HN)=psalic_h_clean$PG.ProteinGroups
DT_limma <- HN %>%
  # Convert NA to 0
  replace(is.na(.), 0) %>%
  mutate(
    # Calculate missing values
    missing_value = rowSums(. == 0),
    missing_value_c = rowSums(select(., matches("HN.*Channel1")) == 0),
    missing_value_t = rowSums(select(., matches("HN.*Channel3")) == 0)
  ) %>%
  filter(
    # Filter rows where missing_value equals the number of columns\
    missing_value != (ncol(select(., matches("HN.*Channel"))))
  )  %>%
  select(-starts_with("missing_value")) %>%  # Remove columns related to missingness counts
  mutate(across(everything(), ~ log2(. + 1)))  # Apply log2 transformation to all columns


group_list <- factor(c(rep('treatment',4),
                       rep("control",4)),
                     levels = c('treatment',"control"))
limma_design <- model.matrix(~0+group_list)
colnames(limma_design) <- levels(group_list)
rownames(limma_design) <- colnames(DT_limma)
cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = limma_design)

#limma
fit <- lmFit(DT_limma, limma_design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, trend=TRUE)

result_limma <- topTable(fit2, coef = 1, n = Inf) %>%
  mutate(
    Group = 'Others',
    Group = if_else(logFC >= 1, 'UP', Group),
    Group = if_else(logFC <= -1, 'DOWN', Group),
    Group = if_else(adj.P.Val >= 0.01, 'Others', Group)
  ) 
result_limma$PG.ProteinGroups <- rownames(result_limma)
merged_data <- psalic_h_clean %>%
  select(1, 2, matches('HN')) %>%
  merge(result_limma, by = "PG.ProteinGroups")
  
write.csv(merged_data[order(merged_data$adj.P.Val),],'pSILAC/DE/Neuri_limma.csv')

ggplot(result_limma, aes(x=AveExpr, y=logFC)) +
  geom_point(aes(color = Group)) +
  scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                     values=c("#67a9cf", "#969696","#ef8a62"),
                     labels = c("Neuri heavy", "not sig", "Neuri light"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_hline(yintercept=-1, linetype="dashed")+ 
  geom_hline(yintercept=1, linetype="dashed")+
  theme_classic()
ggsave('pSILAC/DE/MAplot_Neuri.pdf',width = 6,height = 5)

ggplot(result_limma, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color = Group)) +
  scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                     values=c("#67a9cf", "#969696","#ef8a62"),
                     labels = c("Neuri heavy", "not sig", "Neuri light"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom")+
  geom_hline(yintercept=-log10(0.01), linetype="dashed")+
  geom_vline(xintercept=1, linetype="dashed")+ 
  geom_vline(xintercept=-1, linetype="dashed")+
  theme_classic()
ggsave('pSILAC/DE/vocanol_Neuri.pdf',width = 6,height = 5)

###HS###########

HS=(psalic_h_clean) %>%
  select( matches('HS'))
row.names(HS)=psalic_h_clean$PG.ProteinGroups
DT_limma <- HS %>%
  # Convert NA to 0
  replace(is.na(.), 0) %>%
  mutate(
    # Calculate missing values
    missing_value = rowSums(. == 0),
    missing_value_c = rowSums(select(., matches("HS.*Channel1")) == 0),
    missing_value_t = rowSums(select(., matches("HS.*Channel3")) == 0)
  ) %>%
  filter(
    # Filter rows where missing_value equals the number of columns\
    missing_value != (ncol(select(., matches("HS.*Channel"))))
  )  %>%
  select(-starts_with("missing_value")) %>%  # Remove columns related to missingness counts
  mutate(across(everything(), ~ log2(. + 1)))  # Apply log2 transformation to all columns


group_list <- factor(c(rep('treatment',3),
                       rep("control",3)),
                     levels = c('treatment',"control"))
limma_design <- model.matrix(~0+group_list)
colnames(limma_design) <- levels(group_list)
rownames(limma_design) <- colnames(DT_limma)
cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = limma_design)

#limma
fit <- lmFit(DT_limma, limma_design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, trend=TRUE)

result_limma <- topTable(fit2, coef = 1, n = Inf) %>%
  mutate(
    Group = 'Others',
    Group = if_else(logFC >= 1, 'UP', Group),
    Group = if_else(logFC <= -1, 'DOWN', Group),
    Group = if_else(adj.P.Val >= 0.01, 'Others', Group)
  ) 
result_limma$PG.ProteinGroups <- rownames(result_limma)
merged_data <- psalic_h_clean %>%
  select(1, 2, matches('HS')) %>%
  merge(result_limma, by = "PG.ProteinGroups")

write.csv(merged_data[order(merged_data$adj.P.Val),],'pSILAC/DE/Soma_limma.csv')

ggplot(result_limma, aes(x=AveExpr, y=logFC)) +
  geom_point(aes(color = Group)) +
  scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                     values=c("#67a9cf", "#969696","#ef8a62"),
                     labels = c("Soma heavy", "not sig", "Soma light"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_hline(yintercept=-1, linetype="dashed")+ 
  geom_hline(yintercept=1, linetype="dashed")+
  theme_classic()
ggsave('pSILAC/DE/MAplot_Soma.pdf',width = 6,height = 5)

ggplot(result_limma, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color = Group)) +
  scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                     values=c("#67a9cf", "#969696","#ef8a62"),
                     labels = c("Soma heavy", "not sig", "Soma light"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom")+
  geom_hline(yintercept=-log10(0.01), linetype="dashed")+
  geom_vline(xintercept=1, linetype="dashed")+ 
  geom_vline(xintercept=-1, linetype="dashed")+
  theme_classic()
ggsave('pSILAC/DE/vocanol_Soma.pdf',width = 6,height = 5)




#QuantCAT####
quantcat_H <- quantcat %>%
  select(PG.ProteinGroups, PG.Genes, matches('QH'))%>%
  mutate(across(contains("raw"), ~ ifelse(. < 1000, NA, .)))%>%
  psilac_standardize_format()

quantcat_H_long=melt_intensity_table(quantcat_H)

plot_silac_pg_counts


##data######
QH_L=quantcat_H %>%
  select(1, 2, matches('QH..Channel1'))%>%
  filter(rowSums(is.na(select(., matches("QH..Channel1")))) < (ncol(select(., matches("QH..Channel1"))) / 2))%>%
  rowwise() %>%
  mutate(row_mean = mean(c_across(matches("Channel")), na.rm = TRUE))
QH_H=quantcat_H %>%
  select(1, 2, matches('QH..Channel3'))%>%
  filter(rowSums(is.na(select(., matches("QH..Channel3")))) < (ncol(select(., matches("QH..Channel3"))) / 2))%>%
  rowwise() %>%
  mutate(row_mean = mean(c_across(matches("Channel")), na.rm = TRUE))
QHf_L=quantcat_H %>%
  select(1, 2, matches('QH.f.Channel1'))%>%
  filter(rowSums(is.na(select(., matches("QH.f.Channel1")))) < (ncol(select(., matches("QH.f.Channel1"))) / 2))%>%
  rowwise() %>%
  mutate(row_mean = mean(c_across(matches("Channel")), na.rm = TRUE))
QHf_H=quantcat_H %>%
  select(1, 2, matches('QH.f.Channel3'))%>%
  filter(rowSums(is.na(select(., matches("QH.f.Channel3")))) < (ncol(select(., matches("QH.f.Channel3"))) / 2))%>%
  rowwise() %>%
  mutate(row_mean = mean(c_across(matches("Channel")), na.rm = TRUE))

merged_data_quantcat <- full_join(QH_H, QH_L, by = "Genes", suffix = c(".QH_H", ".QH_L")) %>%
  full_join(QHf_H, by = "Genes") %>%
  full_join(QHf_L, by = "Genes", suffix = c(".QHf_H", ".QHf_L")) %>%
  select(Genes, matches('row_mean'))

write.csv(merged_data_quantcat,'QuantCAT//mean_value.csv')
###plot###########
library(eulerr)
ven_quantcat=list(
  QH_H=unique(na.omit(QH_H$Genes)),
  QH_L=unique(na.omit(QH_L$Genes)),
  QHf_L=unique(na.omit(QHf_L$Genes)),
  QHf_H=unique(na.omit(QHf_H$Genes)))
vd=euler(ven_quantcat)
pdf("QuantCAT//DE/ven_all.pdf",height = 6,width = 7)
par(mar = c(2, 5, 5, 2))
plot(vd, 
     factor_names = TRUE, labels=list(font=2, cex=1.2),
     counts = TRUE,
     key=TRUE,
     cex=1, fills = list(fill = c("#fbb4ae", "#b3cde3", "#ccebc5",'#dfc27d')),
     edges = FALSE,
     quantities = list(type = "counts", "percent"))
dev.off()

ven_quantcat=list(
  QH_H=unique(na.omit(QH_H$Genes)),
  QH_L=unique(na.omit(QH_L$Genes)))
vd=euler(ven_quantcat)
pdf("QuantCAT//DE/QH.pdf",height = 6,width = 7)
par(mar = c(2, 5, 5, 2))
plot(vd, 
     factor_names = TRUE, labels=list(font=2, cex=1.2),
     counts = TRUE,
     key=TRUE,
     cex=1, fills = list(fill = c("#fbb4ae", "#b3cde3")),
     edges = FALSE,
     quantities = list(type = "counts", "percent"))
dev.off()


ven_quantcat=list(
  
  QHf_L=unique(na.omit(QHf_L$Genes)),
  
  QHf_H=unique(na.omit(QHf_H$Genes)))
vd=euler(ven_quantcat)
pdf("QuantCAT//DE/ven_QHf.pdf",height = 6,width = 7)
par(mar = c(2, 5, 5, 2))
plot(vd, 
     factor_names = TRUE, labels=list(font=2, cex=1.2),
     counts = TRUE,
     key=TRUE,
     cex=1, fills = list(fill = c( "#ccebc5",'#dfc27d')),
     edges = FALSE,
     quantities = list(type = "counts", "percent"))
dev.off()


ven_quantcat=list(
  QH_H=unique(na.omit(QH_H$Genes)),
  QHf_H=unique(na.omit(QHf_H$Genes)))
vd=euler(ven_quantcat)
pdf("QuantCAT//DE/QH_QHf_H.pdf",height = 6,width = 7)
par(mar = c(2, 5, 5, 2))
plot(vd, 
     factor_names = TRUE, labels=list(font=2, cex=1.2),
     counts = TRUE,
     key=TRUE,
     cex=1, fills = list(fill = c("#fbb4ae", "#b3cde3")),
     edges = FALSE,
     quantities = list(type = "counts", "percent"))
dev.off()


ven_quantcat=list(
  
  QHf_L=unique(na.omit(QHf_L$Genes)),
  
  QH_L=unique(na.omit(QH_L$Genes)))
vd=euler(ven_quantcat)
pdf("QuantCAT//DE/ven_QHf_QH_L.pdf",height = 6,width = 7)
par(mar = c(2, 5, 5, 2))
plot(vd, 
     factor_names = TRUE, labels=list(font=2, cex=1.2),
     counts = TRUE,
     key=TRUE,
     cex=1, fills = list(fill = c( "#ccebc5",'#dfc27d')),
     edges = FALSE,
     quantities = list(type = "counts", "percent"))
dev.off()
##DE MA and vocanol #######################
###QH##############
QH=data.frame(quantcat_H) %>%
  filter( Genes %in% merged_data_quantcat$Genes)
row.names(QH)=QH$Protein_Group
DT_limma <- QH %>%
  select(matches('QH..Chan'))%>%
  # Convert NA to 0
  replace(is.na(.), 0) %>%
  mutate(
    # Calculate missing values
    missing_value = rowSums(. == 0),
    missing_Channel1 = rowSums(select(., matches("QH..Channel1")) == 0),
    missing_Channel3 = rowSums(select(., matches("QH..Channel3")) == 0)
  ) %>%
  filter(
    # Filter rows where missing_value equals the number of columns\
    missing_value != (ncol(select(., matches("QH..Channel")))),
    missing_Channel1 !=2 ,
    missing_Channel3 !=2,
  )  %>%
  select(-starts_with("missing")) %>%  # Remove columns related to missingness counts
  mutate(across(everything(), ~ log2(. + 1)))  # Apply log2 transformation to all columns


group_list <- factor(c(rep('treatment',3),
                       rep("control",3)),
                     levels = c('treatment',"control"))
limma_design <- model.matrix(~0+group_list)
colnames(limma_design) <- levels(group_list)
rownames(limma_design) <- colnames(DT_limma)
cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = limma_design)

#limma
fit <- lmFit(DT_limma, limma_design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, trend=TRUE)

result_limma <- topTable(fit2, coef = 1, n = Inf) %>%
  mutate(
    Group = 'Others',
    Group = if_else(logFC >= 1, 'UP', Group),
    Group = if_else(logFC <= -1, 'DOWN', Group),
    Group = if_else(adj.P.Val >= 0.01, 'Others', Group)
  ) 
result_limma$Protein_Group <- rownames(result_limma)
result_limma <- quantcat_H %>%
  select(1, 2, matches('QH..Chan')) %>%
  merge(result_limma, by = "Protein_Group")

write.csv(result_limma[order(result_limma$adj.P.Val),],'QuantCAT//DE/QH_L_VS_QH_H_limma.csv')

ggplot(result_limma, aes(x=AveExpr, y=logFC)) +
  geom_point(aes(color = Group)) +
  scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                     values=c("#67a9cf", "#969696","#ef8a62"),
                     labels = c("QH heavy", "not sig", "QH light"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_hline(yintercept=-1, linetype="dashed")+ 
  geom_hline(yintercept=1, linetype="dashed")+
  theme_classic()
ggsave('pSILAC/DE/MAplot_QH.pdf',width = 6,height = 5)

ggplot(result_limma, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color = Group)) +
  scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                     values=c("#67a9cf", "#969696","#ef8a62"),
                     labels = c("QH heavy", "not sig", "QH light"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom")+
  geom_hline(yintercept=-log10(0.01), linetype="dashed")+
  geom_vline(xintercept=1, linetype="dashed")+ 
  geom_vline(xintercept=-1, linetype="dashed")+
  theme_classic()
ggsave('pSILAC/DE/vocanol_QH.pdf',width = 6,height = 5)


###QH vs QHF L##############

DT_limma <- QH %>%
  select(matches('QH..Channel1|QH.f.Channel1'))%>%
  # Convert NA to 0
  replace(is.na(.), 0) %>%
  mutate(
    # Calculate missing values
    missing_value = rowSums(. == 0),
    missing_QH_L = rowSums(select(., matches("QH..Channel1")) == 0),
    missing_QHF_L = rowSums(select(., matches("QH.f.Channel1")) == 0)
  ) %>%
  filter(
    # Filter rows where missing_value equals the number of columns\
    missing_value != (ncol(select(., matches("Channel")))),
    missing_QH_L !=2 ,
    missing_QHF_L !=2,
  )  %>%
  select(-starts_with("missing")) %>%  # Remove columns related to missingness counts
  mutate(across(everything(), ~ log2(. + 1)))  # Apply log2 transformation to all columns


group_list <- factor(c(rep('treatment',3),
                       rep("control",3)),
                     levels = c('treatment',"control"))
limma_design <- model.matrix(~0+group_list)
colnames(limma_design) <- levels(group_list)
rownames(limma_design) <- colnames(DT_limma)
cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = limma_design)

#limma
fit <- lmFit(DT_limma, limma_design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, trend=TRUE)

result_limma <- topTable(fit2, coef = 1, n = Inf) %>%
  mutate(
    Group = 'Others',
    Group = if_else(logFC >= 1, 'UP', Group),
    Group = if_else(logFC <= -1, 'DOWN', Group),
    Group = if_else(adj.P.Val >= 0.01, 'Others', Group)
  ) 
result_limma$Protein_Group <- rownames(result_limma)
result_limma <- data.frame(quantcat_H) %>%
  select(1, 2, matches('QH..Channel1|QH.f.Channel1')) %>%
  merge(result_limma, by = "Protein_Group")

write.csv(result_limma[order(result_limma$adj.P.Val),],'QuantCAT//DE/QH_L_VS_QHf_L_limma.csv')

ggplot(result_limma, aes(x=AveExpr, y=logFC)) +
  geom_point(aes(color = Group)) +
  scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                     values=c("#67a9cf", "#969696","#ef8a62"),
                     labels = c("QH light", "not sig", "QHf light"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_hline(yintercept=-1, linetype="dashed")+ 
  geom_hline(yintercept=1, linetype="dashed")+
  theme_classic()
ggsave('pSILAC/DE/MAplot_QH_L_VS_QHf_L.pdf',width = 6,height = 5)

ggplot(result_limma, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color = Group)) +
  scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                     values=c("#67a9cf", "#969696","#ef8a62"),
                     labels = c("QH light", "not sig", "QHf light"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom")+
  geom_hline(yintercept=-log10(0.01), linetype="dashed")+
  geom_vline(xintercept=1, linetype="dashed")+ 
  geom_vline(xintercept=-1, linetype="dashed")+
  theme_classic()
ggsave('pSILAC/DE/vocanol_QH_L_VS_QHf_L.pdf',width = 6,height = 5)


###QH vs QHF H##############

DT_limma <- QH %>%
  select(matches('QH..Channel3|QH.f.Channel3'))%>%
  # Convert NA to 0
  replace(is.na(.), 0) %>%
  mutate(
    # Calculate missing values
    missing_value = rowSums(. == 0),
    missing_QH_H = rowSums(select(., matches("QH..Channel3")) == 0),
    missing_QHF_H = rowSums(select(., matches("QH.f.Channel3")) == 0)
  ) %>%
  filter(
    # Filter rows where missing_value equals the number of columns\
    missing_value != (ncol(select(., matches("Channel")))),
    missing_QH_H !=2 ,
    missing_QHF_H !=2,
  )  %>%
  select(-starts_with("missing")) %>%  # Remove columns related to missingness counts
  mutate(across(everything(), ~ log2(. + 1)))  # Apply log2 transformation to all columns


group_list <- factor(c(rep('treatment',3),
                       rep("control",3)),
                     levels = c('treatment',"control"))
limma_design <- model.matrix(~0+group_list)
colnames(limma_design) <- levels(group_list)
rownames(limma_design) <- colnames(DT_limma)
cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = limma_design)

#limma
fit <- lmFit(DT_limma, limma_design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, trend=TRUE)

result_limma <- topTable(fit2, coef = 1, n = Inf) %>%
  mutate(
    Group = 'Others',
    Group = if_else(logFC >= 1, 'UP', Group),
    Group = if_else(logFC <= -1, 'DOWN', Group),
    Group = if_else(adj.P.Val >= 0.01, 'Others', Group)
  ) 
result_limma$Protein_Group <- rownames(result_limma)
result_limma <- data.frame(quantcat_H) %>%
  select(1, 2, matches('QH..Channel3|QH.f.Channel3')) %>%
  merge(result_limma, by = "Protein_Group")

write.csv(result_limma[order(result_limma$adj.P.Val),],'QuantCAT//DE/QH_H_VS_QHf_H_limma.csv')

ggplot(result_limma, aes(x=AveExpr, y=logFC)) +
  geom_point(aes(color = Group)) +
  scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                     values=c("#67a9cf", "#969696","#ef8a62"),
                     labels = c("QH heavy", "not sig", "QHf heavy"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_hline(yintercept=-1, linetype="dashed")+ 
  geom_hline(yintercept=1, linetype="dashed")+
  theme_classic()
ggsave('pSILAC/DE/MAplot_QH_H_VS_QHf_H.pdf',width = 6,height = 5)

ggplot(result_limma, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color = Group)) +
  scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                     values=c("#67a9cf", "#969696","#ef8a62"),
                     labels = c("QH heavy", "not sig", "QHf heavy"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom")+
  geom_hline(yintercept=-log10(0.01), linetype="dashed")+
  geom_vline(xintercept=1, linetype="dashed")+ 
  geom_vline(xintercept=-1, linetype="dashed")+
  theme_classic()
ggsave('pSILAC/DE/vocanol_QH_H_VS_QHf_H.pdf',width = 6,height = 5)
