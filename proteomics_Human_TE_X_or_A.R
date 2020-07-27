library(ggplot2)
library(dplyr)
library(ggpubr)

### load expression across 3 layers ###

expression_X_A = readRDS("/Users/agkaessmann/Dropbox/RibosomeProfiling/translatomePaper/Revisions/submission_scripts_and_data/Data/expression_Human_3_layers.rds")

### correlation across layers ###

across_layer_correlation = data.frame(matrix(NA,nrow=6,ncol=3))
names(across_layer_correlation) = c("organ","comparison","rho")
across_layer_correlation$organ = c("Brain","Brain","Liver","Liver","Testis","Testis")
across_layer_correlation$comparison = c("transcriptome2translatome\nvs\ntranscriptome2proteome","transcriptome2translatome\nvs\ntranslatome2proteome")

expression_X_A_Brain = expression_X_A %>% filter(tissue == "Brain")
across_layer_correlation[1,"rho"] = cor.test(expression_X_A_Brain$transcriptome2translatome,expression_X_A_Brain$transcriptome2proteome,method = "spearman")$estimate
across_layer_correlation[2,"rho"] = cor.test(expression_X_A_Brain$transcriptome2translatome,expression_X_A_Brain$translatome2proteome,method = "spearman")$estimate


expression_X_A_Liver = expression_X_A %>% filter(tissue == "Liver")
across_layer_correlation[3,"rho"] =cor.test(expression_X_A_Liver$transcriptome2translatome,expression_X_A_Liver$transcriptome2proteome,method = "spearman")$estimate
across_layer_correlation[4,"rho"] = cor.test(expression_X_A_Liver$transcriptome2translatome,expression_X_A_Liver$translatome2proteome,method = "spearman")$estimate


expression_X_A_Testis = expression_X_A %>% filter(tissue == "Testis")
across_layer_correlation[5,"rho"] = cor.test(expression_X_A_Testis$transcriptome2translatome,expression_X_A_Testis$transcriptome2proteome,method = "spearman")$estimate
across_layer_correlation[6,"rho"] = cor.test(expression_X_A_Testis$transcriptome2translatome,expression_X_A_Testis$translatome2proteome,method = "spearman")$estimate

ggplot(across_layer_correlation, aes(x = comparison, y = rho, fill = organ,col = "black")) +
  geom_bar(stat = "identity",position = "dodge") +
  facet_grid(. ~ organ) +
  scale_fill_manual(values = c("#3399CC", "#339900", "#FF6600" )) +
  scale_color_manual(values = "black") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_rect(fill = "white", size = 2),
    legend.position = "none"
  ) +
  theme(axis.title = element_text(face = "bold", size = 22)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 14)) +
  theme(strip.text.x = element_text(size = 20)) +
  labs(y = expression(rho)) +
  ylim(-0.2,1)


### X vs autosomes ###

expression_X_A_gathered  = expression_X_A %>% gather(key = "from_to", value = "rank_change",c("transcriptome2translatome","transcriptome2proteome","translatome2proteome"))
expression_X_A_gathered$from_to = factor(expression_X_A_gathered$from_to,levels = c("transcriptome2translatome","transcriptome2proteome","translatome2proteome"))
ggplot(expression_X_A_gathered, aes(x = from_to,y = rank_change, fill = tissue, col = X_or_A)) +
  geom_hline(yintercept=0) +
  geom_boxplot(outlier.shape = NA,
               notch=T) +
  facet_grid(. ~ tissue) +
    scale_fill_manual(values = c("#3399CC", "#339900", "#FF6600" )) +
      scale_color_manual(values = c("black" ,"black")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_rect(fill = "white", size = 2)
  ) +
  theme(axis.title = element_text(face = "bold", size = 22)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 14)) +
  theme(strip.text.x = element_text(size = 20)) +
  scale_y_continuous(limits=c(-6000,8000)) +
  stat_compare_means(label = "p.signif",size = 7)

