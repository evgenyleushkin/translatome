library(ggplot2)
library(tidyverse)
library(Mfuzz)

# Spermatocytes, Round spermatids, Elongating spermatids, Spermatozoa
sperm_pallet = c("#FF1400","#FF5A00","#FFA000","#FFE600")

#### set paths to expression file ####
path_to_spermatogenesis_FPKM = "Data/spermatogenesis_FPKM.txt"
CDS_FPKM = read.table(path_to_spermatogenesis_FPKM, header = T)

### tissue-specificity, Tau ###
Tau <- function(x) {
  #  if(any(is.na(x))) stop('Replace NA with 0.')
  #  if(any(x < 0)) stop('Negative values not allowed (log transformed?)')
  x[is.na(x)] = 0
  x[x<0] = 0
  t <- sum(1 - x / max(x)) / (length(x) - 1)
  return(t)
}

max_stage = function(x) {
  return(names(which.max(x)))
}

germ_line = c("spermatocytes", "roundSpermatids", "elongatingSpermatids", "spermatozoa")

####

for(tissue in c("testis",germ_line)) {
  for(data_type in c("TR","RP")) {
    CDS_FPKM[,paste("median",tissue,data_type,sep="_")] = apply(CDS_FPKM[,grepl(tissue,names(CDS_FPKM)) & grepl(data_type,names(CDS_FPKM))],1,median,na.rm = T)
  }
}

CDS_FPKM[,-c(1:3)] = log2(CDS_FPKM[,-c(1:3)])

CDS_FPKM$Tau = apply(CDS_FPKM[,c("median_spermatocytes_RP","median_roundSpermatids_RP","median_elongatingSpermatids_RP","median_spermatozoa_RP")], 1, Tau)
CDS_FPKM$max_stage = as.vector(apply(CDS_FPKM[,c("median_spermatocytes_RP","median_roundSpermatids_RP","median_elongatingSpermatids_RP","median_spermatozoa_RP")],1,max_stage))


# stage-specific genes
Tau_threshold = 0.8
spermatozoa_genes = (CDS_FPKM$Tau>Tau_threshold) & (CDS_FPKM$max_stage == "median_spermatozoa_RP") & CDS_FPKM$median_spermatozoa_RP>1
spermatozoa_genes[is.na(spermatozoa_genes)] = F

roundSpermatids_genes = (CDS_FPKM$Tau>Tau_threshold) & (CDS_FPKM$max_stage == "median_roundSpermatids_RP") & CDS_FPKM$median_roundSpermatids_RP>1
roundSpermatids_genes[is.na(roundSpermatids_genes)] = F

spermatocytes_genes = (CDS_FPKM$Tau>Tau_threshold) & (CDS_FPKM$max_stage == "median_spermatocytes_RP") & CDS_FPKM$median_spermatocytes_RP>1
spermatocytes_genes[is.na(spermatocytes_genes)] = F

elongatingSpermatids_genes = (CDS_FPKM$Tau>Tau_threshold) & (CDS_FPKM$max_stage == "median_elongatingSpermatids_RP") & CDS_FPKM$median_elongatingSpermatids_RP>1
elongatingSpermatids_genes[is.na(elongatingSpermatids_genes)] = F

germ_line_RP_FPKM_names = c("median_spermatocytes_RP","median_roundSpermatids_RP","median_elongatingSpermatids_RP","median_spermatozoa_RP")

somatic_genes = apply(CDS_FPKM[,germ_line_RP_FPKM_names],1,max,na.rm=T) < 0 & CDS_FPKM$median_testis_TR>0
somatic_genes[is.na(somatic_genes)] = F

immature_germ = (somatic_genes == F) & (spermatozoa_genes == F) & (CDS_FPKM$median_spermatozoa_RP<1)
immature_germ[is.na(immature_germ)] = F


CDS_FPKM$cell_type = "rest"
CDS_FPKM[somatic_genes,"cell_type"] = "somatic_genes"
CDS_FPKM[immature_germ,"cell_type"] = "immature_germ"
CDS_FPKM[spermatozoa_genes,"cell_type"] = "spermatozoa_genes"
CDS_FPKM[roundSpermatids_genes,"cell_type"] = "roundSpermatids_genes"
CDS_FPKM[elongatingSpermatids_genes,"cell_type"] = "elongatingSpermatids_genes"
CDS_FPKM[spermatocytes_genes,"cell_type"] = "spermatocytes_genes"



CDS_FPKM$cell_type = factor(CDS_FPKM$cell_type, levels = 
                              c("rest","somatic_genes","immature_germ","spermatocytes_genes","roundSpermatids_genes","elongatingSpermatids_genes","spermatozoa_genes"))


### median expression of specific genes in full testis

median_expressions_specific_genes = matrix(NA,nrow=5,ncol=2,dimnames = list(c("somatic_genes","spermatocytes_genes","roundSpermatids_genes","elongatingSpermatids_genes","spermatozoa_genes"),c("TR","RP")))
for(cell_type in c("somatic_genes","spermatocytes_genes","roundSpermatids_genes","elongatingSpermatids_genes","spermatozoa_genes")) {
  for(data_type in c("TR","RP")) {
    median_expressions_specific_genes[cell_type,data_type] = median(CDS_FPKM[CDS_FPKM$cell_type == cell_type,paste("median","testis",data_type,sep="_")],na.rm=T)
  }
}

## median expression of specific genes in corresponding stages
median_expressions_cell_types = matrix(NA,nrow=4,ncol=2,dimnames = list(germ_line,c("TR","RP")))
for(cell_type in germ_line) {
  for(data_type in c("TR","RP")) {
    median_expressions_cell_types[cell_type,data_type] = median(CDS_FPKM[CDS_FPKM$cell_type == paste(cell_type,"genes",sep="_"),paste("median",cell_type,data_type,sep="_")],na.rm=T)
  }
}

## normalization of gene expression based on stage-specific genes
CDS_FPKM_normalized = CDS_FPKM[CDS_FPKM$chr != "X",c("Gene_ID","CDS_length","chr",names(CDS_FPKM)[grep("median",names(CDS_FPKM))])]
for(cell_type in germ_line) {
  
  CDS_FPKM_normalized[,paste("median",cell_type,"TR",sep="_")] =  CDS_FPKM_normalized[,paste("median",cell_type,"TR",sep="_")] - 
    (median_expressions_cell_types[cell_type,"TR"] - median_expressions_specific_genes[paste(cell_type,"genes",sep="_"),"TR"]) - median_expressions_specific_genes[paste("somatic","genes",sep="_"),"TR"]
  
  CDS_FPKM_normalized[,paste("median",cell_type,"RP",sep="_")] =  CDS_FPKM_normalized[,paste("median",cell_type,"RP",sep="_")] - 
    (median_expressions_cell_types[cell_type,"RP"] - median_expressions_specific_genes[paste(cell_type,"genes",sep="_"),"RP"]) -  median_expressions_specific_genes[paste("somatic","genes",sep="_"),"RP"]  
  
  
}

all_samples = c("testis", germ_line)
for(tissue in all_samples) {
  CDS_FPKM_normalized[,paste("median",tissue,"TE",sep="_")] =   CDS_FPKM_normalized[,paste("median",tissue,"RP",sep="_")] - CDS_FPKM_normalized[,paste("median",tissue,"TR",sep="_")] 
}

#### Mfuzz clustering ####

FPKM_for_ExpressionSet = CDS_FPKM_normalized[,c(15:18)]
row.names(FPKM_for_ExpressionSet) = CDS_FPKM_normalized$Gene_ID

FPKM_for_ExpressionSet[is.na(FPKM_for_ExpressionSet)] = NA
FPKM_for_ExpressionSet[FPKM_for_ExpressionSet == Inf] = NA
FPKM_for_ExpressionSet[FPKM_for_ExpressionSet == -Inf] = NA

minimalSet <- ExpressionSet(assayData=as.matrix(FPKM_for_ExpressionSet))

minimalSet.r = filter.NA(minimalSet, thres = 0)
tmp <- filter.std(minimalSet.r,min.std=0)

minimalSet.s = standardise(minimalSet.r)

mestimate(minimalSet.s) # m = 2.5

cl <- mfuzz(minimalSet.s,c=5,m=2.5)
mfuzz.plot(minimalSet.s,cl=cl,mfrow=c(3,3))

#### plot TE

CDS_FPKM_normalized_TR = cbind(CDS_FPKM_normalized[,c(1,3)],CDS_FPKM_normalized[,grepl("perm",names(CDS_FPKM_normalized)) & 
                                                                                  grepl("median",names(CDS_FPKM_normalized)) &
                                                                                  grepl("_TR",names(CDS_FPKM_normalized))])
CDS_FPKM_normalized_RP = cbind(CDS_FPKM_normalized[,c(1,3)],CDS_FPKM_normalized[,grepl("perm",names(CDS_FPKM_normalized)) & 
                                                                                  grepl("median",names(CDS_FPKM_normalized)) &
                                                                                  grepl("_RP",names(CDS_FPKM_normalized))])

CDS_FPKM_normalized_TE = cbind(CDS_FPKM_normalized[,c(1,3)],CDS_FPKM_normalized[,grepl("perm",names(CDS_FPKM_normalized)) & 
                                                                                  grepl("median",names(CDS_FPKM_normalized)) &
                                                                                  grepl("_TE",names(CDS_FPKM_normalized))])
CDS_FPKM_normalized_TE$layer = "TE"
names(CDS_FPKM_normalized_TE) = c("Gene_ID", "chr", germ_line,"layer")


CDS_FPKM_normalized_TE$cluster = NA
CDS_FPKM_normalized_TE[CDS_FPKM_normalized_TE$Gene_ID %in% row.names(minimalSet.s@assayData$exprs), "cluster"] = cl$cluster



CDS_FPKM_normalized_TR$cluster = NA
CDS_FPKM_normalized_TR[CDS_FPKM_normalized_TE$Gene_ID %in% row.names(minimalSet.s@assayData$exprs), "cluster"] = cl$cluster

CDS_FPKM_normalized_formatted_TR_gathered = gather(CDS_FPKM_normalized_TR, key = "cell type", value = "log2FPKM", 1:4)
CDS_FPKM_normalized_formatted_TR_gathered[is.na(CDS_FPKM_normalized_formatted_TR_gathered)] = 0

CDS_FPKM_normalized_formatted_TR_gathered$`cell type` = factor(CDS_FPKM_normalized_formatted_TR_gathered$`cell type`, levels = (c("spermatocytes","roundSpermatids","elongatingSpermatids","spermatozoa") ))

CDS_FPKM_normalized_RP$cluster = NA
CDS_FPKM_normalized_RP[CDS_FPKM_normalized_TE$Gene_ID %in% row.names(minimalSet.s@assayData$exprs), "cluster"] = cl$cluster

CDS_FPKM_normalized_formatted_RP_gathered = gather(CDS_FPKM_normalized_RP, key = "cell type", value = "log2FPKM", 1:4)
CDS_FPKM_normalized_formatted_RP_gathered$`cell type` = factor(CDS_FPKM_normalized_formatted_RP_gathered$`cell type`, levels = (c("spermatocytes","roundSpermatids","elongatingSpermatids","spermatozoa") ))


CDS_FPKM_normalized_formatted_TE_gathered = gather(CDS_FPKM_normalized_TE, key = "cell type", value = "log2FPKM", 3:6)
CDS_FPKM_normalized_formatted_TE_gathered = CDS_FPKM_normalized_formatted_TE_gathered[CDS_FPKM_normalized_formatted_TR_gathered$log2FPKM>0,]
CDS_FPKM_normalized_formatted_TE_gathered$`cell type` = factor(CDS_FPKM_normalized_formatted_TE_gathered$`cell type`, levels = (c("spermatocytes","roundSpermatids","elongatingSpermatids","spermatozoa") ))



ggplot(CDS_FPKM_normalized_formatted_TE_gathered,aes(x = `cell type`,y=log2FPKM, fill = `cell type`)) +
  geom_hline(yintercept = 0,
             linetype = "solid",
             color = "black") +
    geom_violin(width = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, outlier.shape = NA, show.legend = F) +
  scale_fill_manual(values = rev(sperm_pallet)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_rect(fill = "white", size = 2)
  ) +
  theme(axis.title = element_text(face = "bold", size = 22)) +
  theme(legend.text = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(plot.title = element_text(size = 24,hjust = 0.5,vjust = 0.5,face = "bold")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")) +
  labs(x = "", y = bquote('log'[2]~'TE')) +
  coord_cartesian(ylim = c(-7,7))


