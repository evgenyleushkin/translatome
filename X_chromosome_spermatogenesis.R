library(tidyverse)

organisms_list = c("mouse","chicken")
tissues_list = c("spermatocytes","roundSpermatids","elongatingSpermatids","spermatozoa","testis")
data_type_list = c("rna","ribo")

raw_FPKMs = readRDS("Data/expression_for_X_chromosome_spermatogenesis.rds")

chr_name_Mouse_Chicken = read.table(
  paste(
    "Data/ortholog_pairs/mouse_chicken_1to1.v87.txt",
    sep = ""
  ),
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE
)


raw_FPKMs[raw_FPKMs$organism == "chicken","chrMouse"] = chr_name_Mouse_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "chicken","Gene_ID"],chr_name_Mouse_Chicken$Chicken),"ChrMouse"]
raw_FPKMs[raw_FPKMs$organism == "mouse","chrMouse"] = chr_name_Mouse_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "mouse","Gene_ID"],chr_name_Mouse_Chicken$Mouse),"ChrMouse"]


raw_FPKMs[raw_FPKMs$organism == "mouse","chrChicken"] = chr_name_Mouse_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "mouse","Gene_ID"],chr_name_Mouse_Chicken$Mouse),"ChrChicken"]
raw_FPKMs[raw_FPKMs$organism == "chicken","chrChicken"] = chr_name_Mouse_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "chicken","Gene_ID"],chr_name_Mouse_Chicken$Chicken),"ChrChicken"]

#########

raw_FPKMs[is.na(raw_FPKMs$chrMouse),"chrMouse"] = 0
raw_FPKMs[raw_FPKMs$chrMouse != 0,"chrAX"] = "A"
raw_FPKMs[raw_FPKMs$chrMouse == 0,"chrAX"] = 0
raw_FPKMs[raw_FPKMs$chrMouse == "X","chrAX"] = "X"


raw_FPKMs$CDS_FPKM_A = "NA"
raw_FPKMs[raw_FPKMs$chrAX == "A", "CDS_FPKM_A"]  = raw_FPKMs[raw_FPKMs$chrAX == "A", "CDS_FPKM"] 
raw_FPKMs$CDS_FPKM_A = as.numeric(raw_FPKMs$CDS_FPKM_A)

raw_FPKMs = as.tibble(raw_FPKMs)

median_FPKMs = raw_FPKMs %>% filter(chrAX == "A") %>% group_by(organism,tissue,data_type,replicate) %>% summarise(median_FPKM = median(CDS_FPKM,na.rm=T))
raw_FPKMs = raw_FPKMs %>% group_by(organism,tissue,data_type,replicate) %>% mutate(normalised_FPKM = (CDS_FPKM * mean(median_FPKMs$median_FPKM))/median(CDS_FPKM_A,na.rm=T))


normalised_median_FPKMs = raw_FPKMs %>% filter(chrAX == "A") %>% group_by(organism,tissue,data_type,replicate) %>% summarise(median_FPKM = median(normalised_FPKM,na.rm=T))
normalised_median_FPKMs_X = raw_FPKMs %>% filter(chrAX == "X") %>% group_by(organism,tissue,data_type,replicate) %>% summarise(median_FPKM = median(normalised_FPKM,na.rm=T))

processed_FPKMs = raw_FPKMs %>% group_by(Gene_ID,organism,chr,chrMouse,chrChicken,chrAX,tissue,data_type) %>% summarise(FPKM = mean(normalised_FPKM,na.rm=T))

processed_FPKMs$organism = factor(processed_FPKMs$organism, levels = rev(c("mouse","chicken")))
processed_FPKMs$data_type = factor(processed_FPKMs$data_type, levels = (c("rna","ribo")))



### normalise to chicken FPKM every gene ###


for(tissue in tissues_list) {
  for(data_type in c("rna","ribo")) {
    
    #mouse-chicken
    mouse_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "mouse" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],FPKM)
    chicken_matched_genes = chr_name_Mouse_Chicken[match(pull(processed_FPKMs[processed_FPKMs$organism == "mouse" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,],Gene_ID),chr_name_Mouse_Chicken$Mouse),"Chicken"]
    chicken_FPKM = pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == "testis" & processed_FPKMs$data_type == data_type ,],FPKM)
    chicken_matched_FPKM =  chicken_FPKM[match(chicken_matched_genes, pull(processed_FPKMs[processed_FPKMs$organism == "chicken" & processed_FPKMs$tissue == "testis" & processed_FPKMs$data_type == data_type ,],Gene_ID))]
    processed_FPKMs[processed_FPKMs$organism == "mouse" & processed_FPKMs$tissue == tissue & processed_FPKMs$data_type == data_type,"FPKM_chickened"] = mouse_FPKM/chicken_matched_FPKM
    
    
  }
}



#### plot the X to proto-X ratios ####

median_ratios = data.frame(matrix(NA,ncol=4))
median_ratios = median_ratios[-1,]
names(median_ratios) = c("organism","tissue","data_type","median_ratio")

sampled_median_ratios = data.frame(matrix(NA,ncol=4))
sampled_median_ratios = sampled_median_ratios[-1,]
names(sampled_median_ratios) = c("organism","tissue","data_type","median_ratio")



for(organism in "mouse") {
  for (tissue in tissues_list) {
    for (data_type in data_type_list) {
      pulled_data = pull(processed_FPKMs[processed_FPKMs$organism == organism &
                                           processed_FPKMs$tissue == tissue &
                                           processed_FPKMs$data_type == data_type &
                                           processed_FPKMs$chrMouse == "X",], FPKM_chickened)
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

median_ratios$tissue = factor(median_ratios$tissue, levels = rev(tissues_list))
#row.names(median_ratios) = NA

diff_median_ratios = median_ratios[seq(1,30,by = 3),]
diff_median_ratios2plot = diff_median_ratios %>% filter(data_type == "ribo") %>% select(organism,tissue) 
diff_median_ratios = diff_median_ratios %>% filter(data_type == "ribo") %>% select(median_ratio) - diff_median_ratios %>% filter(data_type == "rna") %>% select(median_ratio)


diff_median_ratios2plot$ratio_diff = diff_median_ratios$median_ratio

sperm_pallet = c("#FF1400","#FF5A00","#FFA000","#FFE600")
data_type.labs = c("Transcriptional\nlayer", "Translational\nlayer")
names(data_type.labs) = c("rna", "ribo")
ggplot(median_ratios[median_ratios$tissue != "testis",], aes(y = median_ratio, x= tissue, fill = data_type, color = tissue)) +
  geom_pointrange(
    stat = "summary",
    fun.ymin = min,
    fun.ymax = max,
    fun.y = median,
    show.legend = F,
    #    fill =  c("white", sperm_pallet[1],"white", sperm_pallet[2],"white", sperm_pallet[3],"white", sperm_pallet[4]),
    fill =  c("white","white","white","white" ,sperm_pallet[1], sperm_pallet[2], sperm_pallet[3], sperm_pallet[4]),
    #    color = myPallet[tissue],
    size = 1.5,
    shape = 22
  ) +
  geom_hline(yintercept = -4,
             linetype = "dotted",
             color = "black") +
  geom_hline(yintercept = -3,
             linetype = "dotted",
             color = "black") +
  geom_hline(yintercept = -2,
             linetype = "dotdash",
             color = "black") +
  geom_hline(yintercept = -1,
             linetype = "solid",
             color = "black") +
  geom_hline(yintercept = 0,
             linetype = "twodash",
             color = "black") +
  geom_hline(yintercept = 1,
             linetype = "dotdash",
             color = "black") +
facet_grid(. ~ data_type, labeller = labeller(data_type = data_type.labs)) +
 scale_color_manual(values = c(sperm_pallet)) +
  labs(x = "", y = bquote('log'[2]~'-ratio')) +
  theme_bw() +
 theme(
    #    panel.border =  element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_rect(fill = "white", size = 2)
  ) +
  theme(axis.title = element_text(face = "bold", size = 22)) +
  theme(legend.text = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) +
  #  theme(axis.text.y = element_text(size = 20, face = c("plain","plain","plain","bold","bold","bold"), color = c("black","black","black",buffering_color,buffering_color,buffering_color))) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 24)) +
  ggtitle("X vs proto-X in germ cells") +
  theme(plot.title = element_text(size = 28, face = "bold")) +
  theme(plot.title = element_text(size = 24)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")) +
  coord_flip(ylim=c(-4.5, 1.5)) 

wilcox_values = data.frame(matrix(NA,ncol=3))
wilcox_values = wilcox_values[-1,]
names(wilcox_values) = c("organism","tissue","value")

for(organism in organisms_list[-2]) {
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


