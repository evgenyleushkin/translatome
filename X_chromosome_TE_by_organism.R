library(tidyverse)

organisms_list = c("human","macaque","mouse","opossum","platypus","chicken")
tissues_list = c("brain","liver","testis")
data_type_list = c("rna","ribo")

# Brain, Liver, Testis
myPallet = c("#3399CC","#339900","#FF6600")
names(myPallet) = c("brain","liver","testis")

####

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

####

raw_FPKMs[raw_FPKMs$organism == "chicken","chrHuman"] = chr_name_Human_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "chicken","Gene_ID"],chr_name_Human_Chicken$Chicken),"ChrHuman"]
raw_FPKMs[raw_FPKMs$organism == "macaque","chrHuman"] = chr_name_Human_Macaque[match(raw_FPKMs[raw_FPKMs$organism == "macaque","Gene_ID"],chr_name_Human_Macaque$Macaque),"ChrHuman"]
raw_FPKMs[raw_FPKMs$organism == "mouse","chrHuman"] = chr_name_Human_Mouse[match(raw_FPKMs[raw_FPKMs$organism == "mouse","Gene_ID"],chr_name_Human_Mouse$Mouse),"ChrHuman"]
raw_FPKMs[raw_FPKMs$organism == "opossum","chrHuman"] = chr_name_Human_Opossum[match(raw_FPKMs[raw_FPKMs$organism == "opossum","Gene_ID"],chr_name_Human_Opossum$Opossum),"ChrHuman"]
raw_FPKMs[raw_FPKMs$organism == "platypus","chrHuman"] = chr_name_Human_Platypus[match(raw_FPKMs[raw_FPKMs$organism == "platypus","Gene_ID"],chr_name_Human_Platypus$Platypus),"ChrHuman"]
raw_FPKMs[raw_FPKMs$organism == "human","chrHuman"] = chr_name_Human_Chicken[match(raw_FPKMs[raw_FPKMs$organism == "human","Gene_ID"],chr_name_Human_Chicken$Human),"ChrHuman"]


raw_FPKMs[is.na(raw_FPKMs$chrHuman),"chrHuman"] = 0
raw_FPKMs[raw_FPKMs$chrHuman != 0,"chrAX"] = "A"
raw_FPKMs[raw_FPKMs$chrHuman == 0,"chrAX"] = 0
raw_FPKMs[raw_FPKMs$chrHuman == "X","chrAX"] = "X"


processed_FPKMs = raw_FPKMs %>% group_by(Gene_ID,organism,tissue,chr,data_type,replicate) %>% 
  filter(is.na(CDS_FPKM) ==  F) %>% 
  filter(CDS_FPKM > 1) %>% 
  mutate(logFPKM = log2(CDS_FPKM)) %>% 
  dplyr::select(-CDS_FPKM) %>% 
  dplyr::select(-CDS_read_count) %>% 
  spread(data_type,logFPKM)  %>% 
  mutate(logTE = ribo - rna)


median_logTE = processed_FPKMs %>% filter(chrAX == "A") %>% group_by(organism,tissue,replicate) %>% dplyr::summarise(median_logTE = median(logTE,na.rm=T))

joined_FPKMs = inner_join(processed_FPKMs,median_logTE)
processed_FPKMs$normalized_TE = processed_FPKMs$logTE - joined_FPKMs$median_logTE

processed_FPKMs$organism = factor(processed_FPKMs$organism, levels = organisms_list)

ggplot(processed_FPKMs %>% filter(chrAX != 0),aes(x = organism,y = normalized_TE,fill = chrAX)) +
  geom_hline(yintercept=0) +
  geom_boxplot(outlier.shape = NA,
               notch=T) +
  facet_grid(. ~ tissue) +
  scale_fill_manual(values = c("grey","red")) +
  theme_bw() +
  theme(
    #    panel.border =  element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_rect(fill = "white", size = 2),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28)
  ) +
  theme(axis.title = element_text(face = "bold", size = 32)) +
  theme(axis.text.y = element_text(size = 24)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 30)) +
  theme(strip.text.x = element_text(size = 30)) +
  scale_y_continuous(limits=c(-4,4)) 


