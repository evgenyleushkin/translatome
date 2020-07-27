library(ggplot2)
library(ape)
library(tidyverse)

organisms_list = c("Human","Macaque","Mouse","Opossum","Platypus","Chicken")
myPallet = c("#3399CC","#339900","#FF6600")

#### set paths to expression file ####

expression_filtered = read.table("Data/expression_filtered.txt",header = T)

#### estimate of within-species vriance and measurement error ####

error_estimate = tibble( colname = colnames(expression_filtered) ) %>%
  filter( ! str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
  group_by( species, tissue, type ) %>%
  summarise( meanVar = mean( genefilter::rowVars( expression_filtered[,colname,drop=FALSE] ), na.rm=TRUE ) ) 
error_estimate = data.frame(error_estimate)

error_estimate$dividedVar = error_estimate$meanVar/2
error_estimate[error_estimate$species == "Human" & error_estimate$tissue == "Testis","dividedVar"] = error_estimate[error_estimate$species == "Human" & error_estimate$tissue == "Testis","meanVar"]
error_estimate[error_estimate$species == "Platypus" & error_estimate$tissue == "Brain" & error_estimate$type == "TR","dividedVar"] = error_estimate[error_estimate$species == "Platypus" & error_estimate$tissue == "Brain" & error_estimate$type == "TR","meanVar"]

####  var estimates ####
var_estimates = tibble( colname = colnames(expression_filtered) ) %>%
  filter( str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
  group_by_all() %>%
  summarise( var = var( expression_filtered[ , colname ], na.rm=TRUE ) ) %>%
  ungroup 
var_estimates = data.frame(var_estimates)

var_estimates$var_corrected = var_estimates$var - error_estimate$dividedVar # correction for within-species vriance and measurement error
var_estimates$species = factor(var_estimates$species, levels = organisms_list)
var_estimates[var_estimates$type == "RP" & var_estimates$tissue != "Testis", "var_corrected"]/var_estimates[var_estimates$type == "TR"  & var_estimates$tissue != "Testis", "var_corrected"]

#### plot dynamic ranges for 2 expression layers across organs and species ####
ggplot(var_estimates, aes(x = tissue, y = var_corrected, alpha = type, col = tissue)) +
  geom_point(size = 5) +
  facet_grid(. ~ species) +
  ylim(0,10) +
  scale_colour_manual(values=myPallet) +
  scale_alpha_manual(values = c(1,0.4)) +
  theme_bw() +
  theme(
    panel.border =  element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_rect(fill = "white", size = 2)
  )  +
  theme(axis.title = element_text(face = "bold", size = 22)) +
  theme(legend.text = element_text(size = 24)) +
  theme(legend.title = element_text(size = 24)) +
  theme(axis.text.y = element_text(size = 24)) +
  theme(axis.text.x = element_text(angle=90,hjust=1,size = 20)) +
  theme(strip.text.x = element_text(size = 24)) +
  theme(plot.title = element_text(size = 26)) +
  labs(x = paste(""), y = paste("var"))



