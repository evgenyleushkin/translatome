library(tidyverse)

myPallet = c("#3399CC","#339900","#FF6600")
names(myPallet) = c("Brain","Liver","Testis")

### load expression data ###

expression_filtered = read.table("Data/expression_filtered.txt",header = T)
k_for_tissues = c(1.234364,1.193857,1.001133)
names(k_for_tissues) = c("Brain","Liver","Testis")

### select tissue ###

tissue = "Liver"

### Z-score calculation for Delta ###

overall_blowing_coefficient = k_for_tissues[tissue]

TR_for_variance = expression_filtered[,grepl("median",names(expression_filtered)) & grepl(tissue,names(expression_filtered)) & grepl("TR",names(expression_filtered)) ]
RP_for_variance = expression_filtered[,grepl("median",names(expression_filtered)) & grepl(tissue,names(expression_filtered)) & grepl("RP",names(expression_filtered)) ]

expression_filtered$TR_var = apply(TR_for_variance,1,var)
expression_filtered$RP_var = apply(RP_for_variance,1,var)

### delta bootstrap sampling for statistical inference ###
for(i in 1:nrow(expression_filtered)) {
  delta_range = c()
  for (human_rep in c("R1")) { # go through all possible combinations of replicates
    for (macaque_rep in c("R1", "R2", "R3")) {
      for (mouse_rep in c("R1", "R2", "R3")) {
        for (opossum_rep in c("R1", "R2", "R3")) {
          for (platypus_rep in c("R1", "R2", "R3")) {
            for (chicken_rep in c("R1", "R2", "R3")) {
              TR_var_sampled = var(
                c(
                  expression_filtered[i, paste("Human_", tissue, ".TR_", human_rep, sep =
                                                 "")],
                  expression_filtered[i, paste("Macaque_", tissue, ".TR_", macaque_rep, sep =
                                                 "")],
                  expression_filtered[i, paste("Mouse_", tissue, ".TR_", mouse_rep, sep =
                                                 "")],
                  expression_filtered[i, paste("Opossum_", tissue, ".TR_", opossum_rep, sep =
                                                 "")],
                  expression_filtered[i, paste("Platypus_", tissue, ".TR_", platypus_rep, sep =
                                                 "")],
                  expression_filtered[i, paste("Chicken_", tissue, ".TR_", chicken_rep, sep =
                                                 "")]
                )
              )
              
              RP_var_sampled = var(
                c(
                  expression_filtered[i, paste("Human_", tissue, ".RP_", human_rep, sep =
                                                 "")],
                  expression_filtered[i, paste("Macaque_", tissue, ".RP_", macaque_rep, sep =
                                                 "")],
                  expression_filtered[i, paste("Mouse_", tissue, ".RP_", mouse_rep, sep =
                                                 "")],
                  expression_filtered[i, paste("Opossum_", tissue, ".RP_", opossum_rep, sep =
                                                 "")],
                  expression_filtered[i, paste("Platypus_", tissue, ".RP_", platypus_rep, sep =
                                                 "")],
                  expression_filtered[i, paste("Chicken_", tissue, ".RP_", chicken_rep, sep =
                                                 "")]
                )
              )
              
              
              delta_sampled = RP_var_sampled - (overall_blowing_coefficient * TR_var_sampled)
              delta_range = c(delta_range, delta_sampled)
            }
          }
        }
      }
    }
  }
   
# calculate delta  
  median_delta = expression_filtered[i, "RP_var"] - (overall_blowing_coefficient * expression_filtered[i, "TR_var"])
  
# calculate delta Z-score
  expression_filtered[i, "z_score"] = median_delta / sd(delta_range)
}



expression_filtered$TR_median = apply(TR_for_variance,1,median)
expression_filtered$RP_median = apply(RP_for_variance,1,median)


expression_filtered[,"delta_color"] = "grey"
expression_filtered[is.na(expression_filtered$z_score > 1.96),"z_score"] = 0
expression_filtered$p_value = 2 * pnorm(-abs(expression_filtered$z_score))
expression_filtered$p_value_corrected = p.adjust(2 * pnorm(-abs(expression_filtered$z_score)),method = "BH") # multiple testing correction
expression_filtered$delta = expression_filtered$RP_var - (expression_filtered$TR_var * overall_blowing_coefficient)
expression_filtered[expression_filtered$p_value_corrected < 0.1 & expression_filtered$delta > 0,"delta_color"] = "#9F00FF"
expression_filtered[expression_filtered$p_value_corrected < 0.1 & expression_filtered$delta < 0,"delta_color"] = "#A6AE00"
expression_filtered$size = 0.5
expression_filtered$alpha = 0.3
expression_filtered[expression_filtered$p_value_corrected < 0.1, "size"] = 0.8
expression_filtered[expression_filtered$p_value_corrected < 0.1, "alpha"] = 0.5


# trim for plotting
expression_filtered[expression_filtered$delta < (-3),"delta"] = -3
expression_filtered[expression_filtered$delta > 3,"delta"] = 3
expression_filtered[expression_filtered$TR_median > 9,"TR_median"] = 9

full_plot = ggplot(expression_filtered, aes(x = TR_median ,
                                            y = delta)) +
  geom_point(col = expression_filtered$delta_color,alpha = expression_filtered$alpha, size = expression_filtered$size ,shape=21, fill = expression_filtered$delta_color) +
  geom_abline(intercept = 0, slope = 0) +
  theme_bw() +
  theme(panel.border =  element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),strip.background = element_rect(fill="white",size=2))  +
  scale_color_gradientn( colours=c("#339900","#339900","gray","gray","gray","#FF6600","#FF6600"), 
                         breaks=c(-12,-3.5,-1.5,0,1.5,3.5,12), limits=c(-7,7) ) +
  coord_cartesian(ylim = c(-3, 3), xlim = c(0,9)) +
  theme(axis.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 24)) +
  theme(legend.title = element_text(size = 24)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(strip.text.x = element_text(size=16)) +
  theme(plot.title = element_text(size = 26)) +
  labs(x = bquote('log'[2]~'(FPKM + 1)') , y = expression(Delta)) 

median_delta = median(expression_filtered$RP_var - (overall_blowing_coefficient * expression_filtered$TR_var),na.rm = T)
lower_quantile = quantile(expression_filtered$delta, p = 0.25)
higher_quantile = quantile(expression_filtered$delta, p = 0.75)
IQR_delta = IQR(expression_filtered$RP_var - (overall_blowing_coefficient * expression_filtered$TR_var),na.rm = T)


density_plot = ggplot(expression_filtered, aes(x = delta)) +
  geom_density(fill = myPallet[tissue],show.legend=F,alpha = 0.8) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = median_delta, linetype='dashed') +
  theme_bw() +
  theme(panel.border =  element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),strip.background = element_rect(fill="white",size=2)) +
  scale_colour_manual(values=myPallet) +
  theme(axis.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 24)) +
  theme(legend.title = element_text(size = 24)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10,angle = 90, hjust = 1)) +
  theme(strip.text.x = element_text(size=16)) +
  theme(plot.title = element_text(size = 26))+
  xlim(-3,3 ) +
  ylim(0,1.5) +
  coord_flip() +
  labs(x = paste(""), y = paste("density")) +
  theme(plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))


plot_grid(full_plot, density_plot, nrow = 1, align = c("h"), axis = c("b"), rel_widths = c(3,1))

median(expression_filtered$RP_var - (overall_blowing_coefficient * expression_filtered$TR_var),na.rm = T)
mean(expression_filtered$RP_var - (overall_blowing_coefficient * expression_filtered$TR_var),na.rm = T)
IQR(expression_filtered$RP_var - (overall_blowing_coefficient * expression_filtered$TR_var),na.rm = T)
var(expression_filtered$RP_var - (overall_blowing_coefficient * expression_filtered$TR_var),na.rm = T)

