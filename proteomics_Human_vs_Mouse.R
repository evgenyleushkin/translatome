library(binr)
library(ggplot2)
library(dplyr)
library(ggpubr)

myPallet = c("#3399CC","#339900","#FF6600")
names(myPallet) = c("brain","liver","testis")

### change of ranks between species ###

expression_all_filtered = readRDS("Data/expression_Human_vs_Mouse_3_layers.rds")


expression_all_filtered$proteome_rank_Human = rank(expression_all_filtered$proteome_Human)
expression_all_filtered$proteome_rank_Mouse = rank(expression_all_filtered$proteome_Mouse)

expression_all_filtered$TR_rank_Human = rank(expression_all_filtered$TR_Human)
expression_all_filtered$TR_rank_Mouse = rank(expression_all_filtered$TR_Mouse)

expression_all_filtered$RP_rank_Human = rank(expression_all_filtered$RP_Human)
expression_all_filtered$RP_rank_Mouse = rank(expression_all_filtered$RP_Mouse)

expression_all_filtered$proteome_rank_change = expression_all_filtered$proteome_rank_Human - expression_all_filtered$proteome_rank_Mouse
expression_all_filtered$TR_rank_change = expression_all_filtered$TR_rank_Human - expression_all_filtered$TR_rank_Mouse
expression_all_filtered$RP_rank_change = expression_all_filtered$RP_rank_Human - expression_all_filtered$RP_rank_Mouse


cor.test(expression_all_filtered$TR_rank_change,expression_all_filtered$proteome_rank_change,method = "spearman")
cor.test(expression_all_filtered$RP_rank_change,expression_all_filtered$proteome_rank_change,method = "spearman")
cor.test(expression_all_filtered$RP_rank_change,expression_all_filtered$TR_rank_change,method = "spearman")

### fast/slow translatome vs transcriptome evolution ###
expression_all_filtered$rank_delta = rank(abs(expression_all_filtered$RP_rank_change)) + rank(-abs(expression_all_filtered$TR_rank_change))
quantiles = quantile(expression_all_filtered$rank_delta , seq(0, 1, by = 0.1))
expression_all_filtered$rank_delta_category  = cut(expression_all_filtered$rank_delta, quantiles, include.lowest = TRUE,labels = c("slow","20%","30%","40%","50%","60%","70%","80%","90%","fast"))

ggplot(expression_all_filtered %>% filter(rank_delta_category == "slow" | rank_delta_category == "fast"), aes(x = rank_delta_category, y = abs(proteome_rank_change),fill = rank_delta_category)) +
  geom_boxplot(outlier.shape = NA,alpha = 1,size=1.2,notch = T) +
  scale_fill_manual(values = c("#A6AE00","#9F00FF")) +
  theme_bw() +
  theme(
    panel.border =  element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_rect(fill = "white", size = 2),
    legend.position = "none"
  ) +
  xlab("translatome vs transcriptome") +
  ylab("absolute proteome rank change") +
  theme(axis.title = element_text(face = "bold",size = 20)) +
  theme(legend.text = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x=element_text(size = 22)) +
  stat_compare_means(label = "p.signif",size = 7) +
  scale_y_continuous(limits=c(0,4100))

### evolution at proteome layer across gene categories ###

protein_rank_change_by_category = data.frame(matrix(NA,ncol=2,nrow=0))
gene_categories_genes_list = readRDS("Data/gene_categories_genes_list.rds")
gene_categories_list = c("ref","BE","TS_Brain","pLI.h","pLI.l","HI.s","HI.i","Old","Young")

for(gene_category in gene_categories_list) {
  rank_change_category = expression_all_filtered[expression_all_filtered$Mouse %in% gene_categories_genes_list[[gene_category]], "RP_rank_change"]
  protein_rank_change_by_category = rbind(protein_rank_change_by_category, data.frame(matrix(c(
    rank_change_category, rep(gene_category, length(rank_change_category))
  ), ncol = 2)))
}

names(protein_rank_change_by_category) = c("rank_change","gene_category")

protein_rank_change_by_category$rank_change = as.numeric(as.vector(protein_rank_change_by_category$rank_change))
protein_rank_change_by_category$gene_category = factor(protein_rank_change_by_category$gene_category, levels = gene_categories_list)

my_comparisons = list(c("BE","TS_Brain"),c("pLI.h","pLI.l"),c("HI.s","HI.i"),c("Old","Young"))

ggplot(protein_rank_change_by_category,aes(x = gene_category,y = abs(rank_change))) +
  geom_boxplot(outlier.shape = NA,alpha = 1,size=1.2, notch = T, fill = myPallet["brain"]) +
  theme_bw() +
  theme(panel.border =  element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  axis.line = element_line(colour = "black"), strip.background = element_rect(fill="white",size=2)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,size = 26)) +
  theme(axis.title = element_text(face = "bold",size = 30)) +
  theme(legend.text = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 24)) +
  theme(axis.ticks.x =  element_blank()) +
  theme(strip.text.x = element_text(size=20)) +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  xlab("gene category") +
  ylab("proteome rank change") +
  stat_compare_means(comparisons = my_comparisons, label.y = 4000, label = "p.signif",size = 7) +
  scale_y_continuous(limits=c(0,4100))

