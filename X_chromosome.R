library(tidyverse)

organisms_list = c("human","macaque","mouse","opossum","platypus","chicken")
tissues_list = c("brain","liver","testis")
data_type_list = c("rna","ribo")

# Brain, Liver, Testis
myPallet = c("#3399CC","#339900","#FF6600")
names(myPallet) = c("brain","liver","testis")
#### load expression file ###
raw_FPKMs = readRDS("Data/expression_full_data.rds")

##### read ortholog pairs #####
chr_name_Human_Chicken = read.table(
  paste(
    "Data/ortholog_pairs/human_chicken_1to1.v87.txt",
    sep = ""
  ),
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE
)

chr_name_Human_Platypus = read.table(
  paste(
    "Data/ortholog_pairs/human_platypus_1to1.v87.txt",
    sep = ""
  ),
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE
)

chr_name_Human_Opossum = read.table(
  paste(
    "Data/ortholog_pairs/human_opossum_1to1.v87.txt",
    sep = ""
  ),
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE
)

chr_name_Human_Mouse = read.table(
  paste(
    "Data/ortholog_pairs/human_mouse_1to1.v87.txt",
    sep = ""
  ),
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE
)

chr_name_Human_Macaque = read.table(
  paste(
    "Data/ortholog_pairs/human_macaque_1to1.v87.txt",
    sep = ""
  ),
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE
)

chr_name_Macaque_Chicken = read.table(
  paste(
    "Data/ortholog_pairs/macaque_chicken_1to1.v87.txt",
    sep = ""
  ),
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE
)

chr_name_Mouse_Chicken = read.table(
  paste(
    "Data/ortholog_pairs/mouse_chicken_1to1.v87.txt",
    sep = ""
  ),
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE
)

chr_name_Opossum_Chicken = read.table(
  paste(
    "Data/ortholog_pairs/opossum_chicken_1to1.v87.txt",
    sep = ""
  ),
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE
)

chr_name_Platypus_Chicken = read.table(
  paste(
    "Data/ortholog_pairs/platypus_chicken_1to1.v87.txt",
    sep = ""
  ),
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE
)
#############################

raw_FPKMs[raw_FPKMs$organism == "chicken","chrHuman"] = chr_name_Human_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "chicken","Gene_ID"],chr_name_Human_Chicken$Chicken),"ChrHuman"]
raw_FPKMs[raw_FPKMs$organism == "macaque","chrHuman"] = chr_name_Human_Macaque[match(raw_FPKMs[raw_FPKMs$organism == "macaque","Gene_ID"],chr_name_Human_Macaque$Macaque),"ChrHuman"]
raw_FPKMs[raw_FPKMs$organism == "mouse","chrHuman"] = chr_name_Human_Mouse[match(raw_FPKMs[raw_FPKMs$organism == "mouse","Gene_ID"],chr_name_Human_Mouse$Mouse),"ChrHuman"]
raw_FPKMs[raw_FPKMs$organism == "opossum","chrHuman"] = chr_name_Human_Opossum[match(raw_FPKMs[raw_FPKMs$organism == "opossum","Gene_ID"],chr_name_Human_Opossum$Opossum),"ChrHuman"]
raw_FPKMs[raw_FPKMs$organism == "platypus","chrHuman"] = chr_name_Human_Platypus[match(raw_FPKMs[raw_FPKMs$organism == "platypus","Gene_ID"],chr_name_Human_Platypus$Platypus),"ChrHuman"]
raw_FPKMs[raw_FPKMs$organism == "human","chrHuman"] = chr_name_Human_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "human","Gene_ID"],chr_name_Human_Chicken$Human),"ChrHuman"]

raw_FPKMs[raw_FPKMs$organism == "human","chrChicken"] = chr_name_Human_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "human","Gene_ID"],chr_name_Human_Chicken$Human),"ChrChicken"]
raw_FPKMs[raw_FPKMs$organism == "macaque","chrChicken"] = chr_name_Macaque_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "macaque","Gene_ID"],chr_name_Macaque_Chicken$Macaque),"ChrChicken"]
raw_FPKMs[raw_FPKMs$organism == "mouse","chrChicken"] = chr_name_Mouse_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "mouse","Gene_ID"],chr_name_Mouse_Chicken$Mouse),"ChrChicken"]
raw_FPKMs[raw_FPKMs$organism == "opossum","chrChicken"] = chr_name_Opossum_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "opossum","Gene_ID"],chr_name_Opossum_Chicken$Opossum),"ChrChicken"]
raw_FPKMs[raw_FPKMs$organism == "platypus","chrChicken"] = chr_name_Platypus_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "platypus","Gene_ID"],chr_name_Platypus_Chicken$Platypus),"ChrChicken"]
raw_FPKMs[raw_FPKMs$organism == "chicken","chrChicken"] = chr_name_Human_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "chicken","Gene_ID"],chr_name_Human_Chicken$Chicken),"ChrChicken"]


raw_FPKMs[is.na(raw_FPKMs$chrHuman),"chrHuman"] = 0
raw_FPKMs[raw_FPKMs$chrHuman != 0,"chrAX"] = "A"
raw_FPKMs[raw_FPKMs$chrHuman == 0,"chrAX"] = 0
raw_FPKMs[raw_FPKMs$chrHuman == "X","chrAX"] = "X"


raw_FPKMs$CDS_FPKM_A = "NA"
raw_FPKMs[raw_FPKMs$chrAX == "A", "CDS_FPKM_A"]  = raw_FPKMs[raw_FPKMs$chrAX == "A", "CDS_FPKM"] 
raw_FPKMs$CDS_FPKM_A = as.numeric(raw_FPKMs$CDS_FPKM_A)

raw_FPKMs = as.tibble(raw_FPKMs)


median_FPKMs = raw_FPKMs %>% filter(chrAX == "A") %>% group_by(organism,tissue,data_type,replicate) %>% summarise(median_FPKM = median(CDS_FPKM,na.rm=T))
raw_FPKMs = raw_FPKMs %>% group_by(organism,tissue,data_type,replicate) %>% mutate(normalised_FPKM = (CDS_FPKM * mean(median_FPKMs$median_FPKM))/median(CDS_FPKM_A,na.rm=T))

normalised_median_FPKMs = raw_FPKMs %>% filter(chrAX == "A") %>% group_by(organism,tissue,data_type,replicate) %>% summarise(median_FPKM = median(normalised_FPKM,na.rm=T))
normalised_median_FPKMs_X = raw_FPKMs %>% filter(chrAX == "X") %>% group_by(organism,tissue,data_type,replicate) %>% summarise(median_FPKM = median(normalised_FPKM,na.rm=T))

processed_FPKMs = raw_FPKMs %>% group_by(Gene_ID,organism,chr,chrHuman,chrChicken,chrAX,tissue,data_type) %>% summarise(FPKM = mean(normalised_FPKM,na.rm=T))

processed_FPKMs$organism = factor(processed_FPKMs$organism, levels = rev(c("human","macaque","mouse","opossum","platypus","chicken")))
processed_FPKMs$data_type = factor(processed_FPKMs$data_type, levels = (c("rna","ribo")))

### normalise to chicken FPKM every gene ###

for(tissue in c("brain","liver","testis")) {
  for(data_type in c("rna","ribo")) {
    #human-chicken
    human_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "human" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],FPKM)
    chicken_matched_genes = chr_name_Human_Chicken[match(pull(processed_FPKMs[processed_FPKMs$organism == "human" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],Gene_ID),chr_name_Human_Chicken$Human),"Chicken"]
    chicken_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type ,],FPKM)
    chicken_matched_FPKM =  chicken_FPKM[match(chicken_matched_genes, pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type ,],Gene_ID))]
    processed_FPKMs[processed_FPKMs$organism == "human" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,"FPKM_chickened"] = human_FPKM/chicken_matched_FPKM
    
    #macaque-chicken
    macaque_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "macaque" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],FPKM)
    chicken_matched_genes = chr_name_Macaque_Chicken[match(pull(processed_FPKMs[processed_FPKMs$organism == "macaque" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],Gene_ID),chr_name_Macaque_Chicken$Macaque),"Chicken"]
    chicken_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type ,],FPKM)
    chicken_matched_FPKM =  chicken_FPKM[match(chicken_matched_genes, pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type ,],Gene_ID))]
    processed_FPKMs[processed_FPKMs$organism == "macaque" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,"FPKM_chickened"] = macaque_FPKM/chicken_matched_FPKM
    
    #mouse-chicken
    mouse_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "mouse" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],FPKM)
    chicken_matched_genes = chr_name_Mouse_Chicken[match(pull(processed_FPKMs[processed_FPKMs$organism == "mouse" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],Gene_ID),chr_name_Mouse_Chicken$Mouse),"Chicken"]
    chicken_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type ,],FPKM)
    chicken_matched_FPKM =  chicken_FPKM[match(chicken_matched_genes, pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type ,],Gene_ID))]
    processed_FPKMs[processed_FPKMs$organism == "mouse" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,"FPKM_chickened"] = mouse_FPKM/chicken_matched_FPKM
    
    #opossum-chicken
    opossum_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "opossum" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],FPKM)
    chicken_matched_genes = chr_name_Opossum_Chicken[match(pull(processed_FPKMs[processed_FPKMs$organism == "opossum" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],Gene_ID),chr_name_Opossum_Chicken$Opossum),"Chicken"]
    chicken_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type ,],FPKM)
    chicken_matched_FPKM =  chicken_FPKM[match(chicken_matched_genes, pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type ,],Gene_ID))]
    processed_FPKMs[processed_FPKMs$organism == "opossum" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,"FPKM_chickened"] = opossum_FPKM/chicken_matched_FPKM
    
    #platypus-chicken
    platypus_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "platypus" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],FPKM)
    chicken_matched_genes = chr_name_Platypus_Chicken[match(pull(processed_FPKMs[processed_FPKMs$organism == "platypus" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],Gene_ID),chr_name_Platypus_Chicken$Platypus),"Chicken"]
    chicken_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type ,],FPKM)
    chicken_matched_FPKM =  chicken_FPKM[match(chicken_matched_genes, pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type ,],Gene_ID))]
    processed_FPKMs[processed_FPKMs$organism == "platypus" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,"FPKM_chickened"] = platypus_FPKM/chicken_matched_FPKM
    
  }
}


#### plot the X to proto-X ratios

median_ratios = data.frame(matrix(NA,ncol=4))
median_ratios = median_ratios[-1,]
names(median_ratios) = c("organism","tissue","data_type","median_ratio")

sampled_median_ratios = data.frame(matrix(NA,ncol=4))
sampled_median_ratios = sampled_median_ratios[-1,]
names(sampled_median_ratios) = c("organism","tissue","data_type","median_ratio")

for(organism in organisms_list[-6]) {
  for (tissue in tissues_list) {
    for (data_type in data_type_list) {
      pulled_data = dplyr::pull(processed_FPKMs[processed_FPKMs$organism == organism &
                                                  processed_FPKMs$tissue == tissue &
                                                  processed_FPKMs$data_type == data_type &
                                                  processed_FPKMs$chrHuman == "X",], FPKM_chickened)
      median_ratio = median(pulled_data, na.rm = T)
      bootstrap_replicates = 100
      sampled_data = c()
      for (i in 1:bootstrap_replicates) {
        sampled_data = c(sampled_data, median(sample(
          pulled_data,
          size = length(pulled_data),
          replace = T
        ), na.rm = T))
      }
      
      data_to_add = data.frame(rep(organism, 3),
                               rep(tissue, 3),
                               rep(data_type, 3),
                               c(log2(median_ratio), log2(quantile(sampled_data,0.025)),log2(quantile(sampled_data,0.975))))
      names(data_to_add)  = c("organism","tissue","data_type","median_ratio")
      median_ratios = rbind(median_ratios, data_to_add)
      
      data_to_add = data.frame(rep(organism, bootstrap_replicates),
                               rep(tissue, bootstrap_replicates),
                               rep(data_type, bootstrap_replicates),
                               log(sampled_data))
      names(data_to_add)  = c("organism","tissue","data_type","median_ratio")
      sampled_median_ratios = rbind(sampled_median_ratios, data_to_add)
      
    }
  }
}

###
library(plyr)
###

median_ratios$organism = factor(median_ratios$organism, levels = rev(organisms_list[-6]))
median_ratios$organism  = revalue(median_ratios$organism, c("platypus" = "Platypus","opossum" = "Opossum","mouse" = "Mouse", "macaque" = "Macaque", "human" = "Human"))



tissue = "testis" # select tissue to plot

diff_median_ratios = median_ratios[seq(1,90,by = 3),]
diff_median_ratios2plot = diff_median_ratios %>% filter(data_type == "ribo") %>% select(organism,tissue) 
diff_median_ratios = diff_median_ratios %>% filter(data_type == "ribo") %>% select(median_ratio) - diff_median_ratios %>% filter(data_type == "rna") %>% select(median_ratio)


diff_median_ratios2plot$ratio_diff = diff_median_ratios$median_ratio

data_type.labs = c("Transcriptional\nlayer", "Translational\nlayer")
names(data_type.labs) = c("rna", "ribo")
ggplot(median_ratios[median_ratios$tissue == tissue,], aes(y = median_ratio, x= organism, fill = data_type)) +
  geom_pointrange(
    stat = "summary",
    fun.ymin = min,
    fun.ymax = max,
    fun.y = median,
    show.legend = F,
    color = myPallet[tissue],
    size = 1.5,
    shape = 22
  ) +
  geom_hline(yintercept = -1,
             linetype = "dashed",
             color = "black") +
  geom_hline(yintercept = 0,
             linetype = "solid",
             color = "black") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  geom_vline(xintercept = 1,
             linetype = "dotted",
             color = myPallet[tissue]) +
  geom_vline(xintercept = 2,
             linetype = "dotted",
             color = myPallet[tissue]) +
  geom_vline(xintercept = 3,
             linetype = "dotted",
             color = myPallet[tissue]) +
  geom_vline(xintercept = 4,
             linetype = "dotted",
             color = myPallet[tissue]) +
  geom_vline(xintercept = 5,
             linetype = "dotted",
             color = myPallet[tissue]) +
  facet_grid(. ~ data_type, labeller = labeller(data_type = data_type.labs)) +
  scale_fill_manual(values = c("white", unname(myPallet[tissue]))) +
  labs(x = "", y = bquote('log'[2]~'-ratio')) +
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
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 24)) +
  ggtitle(paste("X vs proto-X in",tolower(tissue))) +
  theme(plot.title = element_text(size = 28, face = "bold")) +
  theme(plot.title = element_text(size = 24)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(plot.margin = unit(c(0.1, 3.5, 0.1, 0.1), "cm")) +
  coord_flip(ylim=c(-2, 2),clip = 'off') 

### wilcoxon test

wilcox_values = data.frame(matrix(NA,ncol=3))
wilcox_values = wilcox_values[-1,]
names(wilcox_values) = c("organism","tissue","value")

for(organism in organisms_list[-6]) {
  for (tissue in tissues_list) {
    wilcox_value  = wilcox.test(sampled_median_ratios[sampled_median_ratios$organism == organism & sampled_median_ratios$tissue == tissue & sampled_median_ratios$data_type == "rna","median_ratio"],
                                sampled_median_ratios[sampled_median_ratios$organism == organism & sampled_median_ratios$tissue == tissue & sampled_median_ratios$data_type == "ribo","median_ratio"])$p.value
    
    data_to_add = data.frame(rep(organism, 1),
                             rep(tissue, 1),
                             wilcox_value)
    names(data_to_add)  = c("organism","tissue","value")
    wilcox_values = rbind(wilcox_values, data_to_add)
  }
}


