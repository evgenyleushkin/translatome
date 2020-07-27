library(tidyverse)
library(ggplot2)
library(scales)

#### load expression data ####
expression_filtered = read.table("Data/expression_filtered.txt",header = T)

####  estimate of within-species variance and measurment error variance ####

error_estimate = tibble( colname = colnames(expression_filtered) ) %>%
  filter( ! str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
  group_by( species, tissue, type ) %>%
  summarise( meanVar = mean( genefilter::rowVars( expression_filtered[,colname,drop=FALSE] ), na.rm=TRUE ) ) 
error_estimate = data.frame(error_estimate)

error_estimate$dividedVar = error_estimate$meanVar/3
error_estimate[error_estimate$species == "Human" & error_estimate$tissue == "Testis","dividedVar"] = error_estimate[error_estimate$species == "Human" & error_estimate$tissue == "Testis","meanVar"]/2
error_estimate[error_estimate$species == "Platypus" & error_estimate$tissue == "Brain" & error_estimate$type == "TR","dividedVar"] = error_estimate[error_estimate$species == "Platypus" & error_estimate$tissue == "Brain" & error_estimate$type == "TR","meanVar"]/2


#### set tissue: Brain, Liver or Testis ####
tissue = "Brain" 
organisms_list = c("Human","Macaque","Mouse","Opossum","Platypus","Chicken")
if(tissue == "Liver") {
  organisms_list = c("Macaque","Mouse","Opossum","Platypus","Chicken")
}


####  variance estimates #####

variance_estimates_RP_TR_TE = tibble( colname = colnames(expression_filtered) ) %>%
  filter( str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
  spread( type, colname ) %>%
  group_by_all() %>%
  summarise( RP_var = var( expression_filtered[ , RP ], na.rm=TRUE ) ) %>%
  ungroup %>% 
  group_by_all() %>%
  summarise( TR_var = var( expression_filtered[ , TR ], na.rm=TRUE ) ) %>%
  ungroup %>%
  group_by_all() %>%
  summarise( TE_var = var( expression_filtered[ , RP ] - expression_filtered[ , TR ], na.rm=TRUE ) ) %>%
  ungroup 
variance_estimates = data.frame(variance_estimates_RP_TR_TE)

#### parameter estimation based on Mouse and Macaque Brain ####
## estimate of standard deviation of the ancestral transcript abundance ##
TR_sd = as.vector(unlist(sqrt((variance_estimates %>% filter(tissue == "Brain", species == "Macaque") %>% select(TR_var) +
                                 variance_estimates %>% filter(tissue == "Brain", species == "Mouse") %>% select(TR_var) - 
                                 (error_estimate %>% filter(tissue == "Brain", species == "Macaque", type == "TR") %>% select(dividedVar) +
                                    error_estimate %>% filter(tissue == "Brain", species == "Mouse", type == "TR") %>% select(dividedVar)))/2)))

## estimate of within-species variance and measurment error variance, transcriptome layer ##
TR_error_sd = as.vector(unlist(sqrt((error_estimate %>% filter(tissue == "Brain", species == "Macaque", type == "TR") %>% select(dividedVar) +
                                       error_estimate %>% filter(tissue == "Brain", species == "Mouse", type == "TR") %>% select(dividedVar))/2)))

## estimate of within-species variance and measurment error variance, translatome layer ##
RP_error_sd = as.vector(unlist(sqrt((error_estimate %>% filter(tissue == "Brain", species == "Macaque", type == "RP") %>% select(dividedVar) +
                                       error_estimate %>% filter(tissue == "Brain", species == "Mouse", type == "RP") %>% select(dividedVar))/2)))

## estimate of standard deviation of the ancestral tranlsation efficiency ##
TE_sd = as.vector(unlist(sqrt((variance_estimates %>% filter(tissue == "Brain", species == "Macaque") %>% select(TE_var) +
                                 variance_estimates %>% filter(tissue == "Brain", species == "Mouse") %>% select(TE_var) - 
                                 (error_estimate %>% filter(tissue == "Brain", species == "Macaque", type == "TR") %>% select(dividedVar) +
                                    error_estimate %>% filter(tissue == "Brain", species == "Mouse", type == "TR") %>% select(dividedVar))-
                                 (error_estimate %>% filter(tissue == "Brain", species == "Macaque", type == "RP") %>% select(dividedVar) +
                                    error_estimate %>% filter(tissue == "Brain", species == "Mouse", type == "RP") %>% select(dividedVar)))/2)))

## estimate of standard deviation of transcriptome layer expression change ##
Macaque_Mouse_TR_change_sd_estimate = sqrt(var((expression_filtered$median_Mouse_Brain.TR - expression_filtered$median_Macaque_Brain.TR),na.rm=T)/2 - (error_estimate %>% filter(tissue == "Brain", species == "Macaque", type == "TR") %>% select(dividedVar) +
                                                                                                                                                         error_estimate %>% filter(tissue == "Brain", species == "Mouse", type == "TR") %>% select(dividedVar))/2 )

## estimate of standard deviation of translatome layer expression change ##
Macaque_Mouse_RP_var_div = (var((expression_filtered$median_Mouse_Brain.RP - expression_filtered$median_Macaque_Brain.RP),na.rm=T) - (error_estimate %>% filter(tissue == "Brain", species == "Macaque", type == "TR") %>% select(dividedVar) +
                                                                                                                                                         error_estimate %>% filter(tissue == "Brain", species == "Mouse", type == "TR") %>% select(dividedVar)))

## calculate expression divergence between macaque and mouse at the translatome layer ##
Macaque_Mouse_RP_var_denumerator = (var(expression_filtered$median_Mouse_Brain.RP) + var(expression_filtered$median_Macaque_Brain.RP))/2 - (error_estimate %>% filter(tissue == "Brain", species == "Macaque", type == "TR") %>% select(dividedVar) +
                                                                                                                                              error_estimate %>% filter(tissue == "Brain", species == "Mouse", type == "TR") %>% select(dividedVar))

RP_real_divergence = as.vector(unlist(Macaque_Mouse_RP_var_div/Macaque_Mouse_RP_var_denumerator))

#### simulate the data ####

## ancestral transcript abundance 
TR_real = rnorm(5060, sd = TR_sd) 
## ancestral TE 
TE_real = rnorm(5060, sd = TE_sd) 

## TE change parameter
TE_diff_range = c(1:60 / 100) 
## compensatory evolution parameter
compensation_range = c((-(20:1) / 100), 0, (1:70) / 100) 
## transcript abundance change
TR_change_sd = Macaque_Mouse_TR_change_sd_estimate[1,1] 
TR_change1 = rnorm(5060, sd = TR_change_sd) 
TR_change2 = rnorm(5060, sd = TR_change_sd) 

## modern transcript abundance
TR_no_errors_S1 = TR_real + TR_change1
TR_no_errors_S2 = TR_real + TR_change2

TR_sim_S1 = TR_no_errors_S1 + rnorm(5060, sd = TR_error_sd)
TR_sim_S2 = TR_no_errors_S2 + rnorm(5060, sd = TR_error_sd)

## divergence estimate
TR_distance = var(TR_sim_S2 - TR_sim_S1)
TR_distance_corrected = matrix(NA, nrow = length(TE_diff_range), ncol = length(compensation_range))
RP_distance = matrix(NA, nrow = length(TE_diff_range), ncol = length(compensation_range))
RP_distance_corrected = matrix(NA, nrow = length(TE_diff_range), ncol = length(compensation_range))
real_distance = matrix(NA, nrow = length(TE_diff_range), ncol = length(compensation_range))


RP_distance_simulated = data.frame(matrix(NA,nrow=0,ncol=3))
names(RP_distance_simulated) = c("Compensation","TE_change","Divergence")

row_num = 1
for(TE_diff in 1:length(TE_diff_range)) { ## vary trhough TE change parameter
  
  TE_change_sd = sqrt(TE_sd^2 * TE_diff_range[TE_diff])
  
  for (compensation in 1:length(compensation_range)) { ## vary trhough compensatory evolution parameter
    
    TE_change1 = rnorm(5060, sd = TE_change_sd)  - (compensation_range[compensation] * TR_change1)
    TE_change2 = rnorm(5060, sd = TE_change_sd)  - (compensation_range[compensation] * TR_change2)
    
    RP_sim_S1 = TR_no_errors_S1 + TE_real + TE_change1 + rnorm(5060, sd = RP_error_sd)
    RP_sim_S2 = TR_no_errors_S2 + TE_real + TE_change2 + rnorm(5060, sd = RP_error_sd)
    
    RP_distance[TE_diff,compensation] = var(RP_sim_S2 - RP_sim_S1)
    RP_distance_corrected[TE_diff,compensation] = (var(RP_sim_S2 - RP_sim_S1) - (RP_error_sd^2))/((var(RP_sim_S1) - RP_error_sd^2 + var(RP_sim_S2) - RP_error_sd^2)/2)
    TR_distance_corrected[TE_diff,compensation] = (var(TR_sim_S2 - TR_sim_S1) - (TR_error_sd^2))/((var(TR_sim_S1) - TR_error_sd^2 + var(TR_sim_S2) - TR_error_sd^2)/2)
    real_distance[TE_diff,compensation] = RP_real_divergence
    
    RP_distance_simulated[row_num,"Compensation"] = compensation_range[compensation]
    RP_distance_simulated[row_num,"TE_change"] = TE_diff_range[TE_diff]
    RP_distance_simulated[row_num,"Divergence"]  = (var(RP_sim_S2 - RP_sim_S1) - (RP_error_sd^2))/((var(RP_sim_S1) - RP_error_sd^2 + var(RP_sim_S2) - RP_error_sd^2)/2)
    
    row_num = row_num + 1 
  }
}










