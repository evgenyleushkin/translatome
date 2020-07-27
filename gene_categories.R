require(methods)
library(ape)
library(tidyverse)

on_desktop = "/Users/agkaessmann/bwFor"
organisms_list = c("Macaque","Mouse","Opossum")

myPallet = c("#3399CC","#339900","#FF6600")
names(myPallet) = c("Brain", "Liver", "Testis")

#### select tissue of interest ####
tissue = "Brain"

#### load expression data ####
all_expression_rmo = readRDS("Data/all_expression_rmo.rds")

#### load gene categories list ####
gene_categories_genes_list = readRDS("Data/gene_categories_genes_list.rds")
gene_categories = c(names(gene_categories_genes_list)[1:6],paste("TS",tissue,sep="_"),names(gene_categories_genes_list)[10:11])

#### filter gene expression ####
median_expression_overall_rmo = apply(all_expression_rmo[,-1],1,median)
expression_filtered_overall = all_expression_rmo[median_expression_overall_rmo>1,] # filter by median expression > 1
expression_filtered_overall[,-1] = log(expression_filtered_overall[,-1] + 1)

####  estimate of within-species variance and measurment error variance ####
error_estimate = tibble( colname = colnames(expression_filtered_overall) ) %>%
  filter( ! str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
  group_by( species, tissue, type ) %>%
  summarise( meanVar = mean( genefilter::rowVars( expression_filtered_overall[,colname,drop=FALSE] ), na.rm=TRUE ) ) 
error_estimate = data.frame(error_estimate)
error_estimate$dividedVar = error_estimate$meanVar/3 ## divide by the number of replicates

####  variance estimates #####
variance_estimates = tibble( colname = colnames(expression_filtered_overall) ) %>%
  filter( str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
  #  filter( tissue == tissue) %>%
  spread( type, colname ) %>%
  group_by_all() %>%
  summarise( RP_var = var( expression_filtered_overall[ , RP ], na.rm=TRUE ) ) %>%
  ungroup %>% 
  group_by_all() %>%
  summarise( TR_var = var( expression_filtered_overall[ , TR ], na.rm=TRUE ) ) %>%
  ungroup %>%
  group_by_all() %>%
  summarise( TE_var = var( expression_filtered_overall[ , RP ] - expression_filtered_overall[ , TR ], na.rm=TRUE ) ) %>%
  ungroup 

variance_estimates = data.frame(variance_estimates)

#### tree lentgth estimate by category and bootstrapping ####
bootstrap.reps = 1000

TR_lengths_bootstrapped = data.frame(matrix(NA, nrow = 1001, ncol = 9))
RP_lengths_bootstrapped = data.frame(matrix(NA, nrow = 1001, ncol = 9))
names(TR_lengths_bootstrapped) = gene_categories
names(RP_lengths_bootstrapped) = gene_categories

for(gene_category in gene_categories) {
  
  data = expression_filtered_overall[expression_filtered_overall$Gene_ID %in% gene_categories_genes_list[[gene_category]],] # subset for genes in a category
  if(gene_category == paste("TS",tissue,sep="_")) {
    median_expression_tissue_rmo = apply(all_expression_rmo[,grepl(tissue,names(all_expression_rmo))],1,median)
    
    expression_filtered_tissue = all_expression_rmo[median_expression_tissue_rmo>1,] # filter by median expression > 1 in a tissue of interest (for tissue-specific genes)
    expression_filtered_tissue[,-1] = log(expression_filtered_tissue[,-1])
    expression_filtered_tissue[expression_filtered_tissue == -Inf] = NA
    
    data = expression_filtered_tissue[expression_filtered_tissue$Gene_ID %in% gene_categories_genes_list[[gene_category]],]
  }
  
  ## proceed first with the actual estimate within category ##
  actual.data = data
  
  ##  estimate of within-species variance and measurment error variance, within category ##
  error_estimate_category = tibble( colname = colnames(actual.data) ) %>%
    filter( ! str_detect( colname, "median" ) ) %>%
    filter( colname != "Gene_ID" ) %>%
    separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
    group_by( species, tissue, type ) %>%
    summarise( meanVar = mean( genefilter::rowVars( actual.data[,colname,drop=FALSE] ), na.rm=TRUE ) ) 
  error_estimate_category = data.frame(error_estimate_category)
  error_estimate_category$dividedVar = error_estimate_category$meanVar/3
  
  ## transcriptome distance measurment, within category ##
  TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>% ## correct for within-species variance and measurment error variance
    filter( sp1 != sp2 ) %>%
    rowwise() %>%
    group_by_all() %>%
    summarise(var_divergence = 
                var(actual.data[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] - 
                      actual.data[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) - 
                error_estimate_category[error_estimate_category$species == sp1 & error_estimate_category$tissue == tissue & error_estimate_category$type == "TR", "dividedVar"] -
                error_estimate_category[error_estimate_category$species == sp2 & error_estimate_category$tissue == tissue & error_estimate_category$type == "TR", "dividedVar"]) %>%
    ungroup %>%
    spread( sp2, var_divergence, fill=0 ) %>%
    column_to_rownames( "sp1" ) %>%
    as.matrix
  
  TR_distances_with_errors_normalized = TR_distances_with_errors ## normalize to dynamic range on transcriptome layer
  for(sp1 in organisms_list) {
    for(sp2 in organisms_list) {
      denominator = variance_estimates[variance_estimates$species == sp1 & variance_estimates$tissue == tissue, "TR_var"] + 
        variance_estimates[variance_estimates$species == sp2 & variance_estimates$tissue == tissue, "TR_var"] -
        error_estimate[error_estimate$species == sp1 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "dividedVar"] -
        error_estimate[error_estimate$species == sp2 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "dividedVar"] 
      TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/(denominator/2)
    }
  }
  
  TR_tree = nj( TR_distances_with_errors_normalized )
  TR_lengths_bootstrapped[1,gene_category] = sum(TR_tree$edge.length)
  
  ## translatome distance measurment, within category ##
  RP_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>% ## correct for within-species variance and measurment error variance
    filter( sp1 != sp2 ) %>%
    rowwise() %>%
    group_by_all() %>%
    summarise(var_divergence = 
                var(actual.data[ , paste("median_",sp1,"_",tissue,".","RP",sep = "")] - 
                      actual.data[ , paste("median_",sp2,"_",tissue,".","RP",sep = "")], na.rm=TRUE ) - 
                error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "RP", "dividedVar"] -
                error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "RP", "dividedVar"]) %>%
    ungroup %>%
    spread( sp2, var_divergence, fill=0 ) %>%
    column_to_rownames( "sp1" ) %>%
    as.matrix
  
  RP_distances_with_errors_normalized = RP_distances_with_errors ## normalize to dynamic range on translatome layer
  for(sp1 in organisms_list) {
    for(sp2 in organisms_list) {
      denominator = variance_estimates[variance_estimates$species == sp1 & variance_estimates$tissue == tissue, "RP_var"] + 
        variance_estimates[variance_estimates$species == sp2 & variance_estimates$tissue == tissue, "RP_var"] -
        error_estimate[error_estimate$species == sp1 & error_estimate$type == "RP" & error_estimate$tissue == tissue, "dividedVar"] -
        error_estimate[error_estimate$species == sp2 & error_estimate$type == "RP" & error_estimate$tissue == tissue, "dividedVar"] 
      RP_distances_with_errors_normalized[sp1,sp2] = RP_distances_with_errors[sp1,sp2]/(denominator/2)
    }
  }
  
  RP_tree = nj( RP_distances_with_errors_normalized )
  RP_lengths_bootstrapped[1,gene_category] = sum(RP_tree$edge.length)
  
  ## bootstrapping, sample with replacament category data ##
  
  for(i in 2:(bootstrap.reps+1)) {
    
    sampled.data <- data[sample(x = 1:dim(data)[1], size = dim(data)[1], replace = TRUE), ]
    
    ##  estimate of within-species variance and measurment error variance, within category, sampled ##
    error_estimate_sampled = tibble( colname = colnames(sampled.data) ) %>%
      filter( ! str_detect( colname, "median" ) ) %>%
      filter( colname != "Gene_ID" ) %>%
      separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
      group_by( species, tissue, type ) %>%
      summarise( meanVar = mean( genefilter::rowVars( sampled.data[,colname,drop=FALSE] ), na.rm=TRUE ) ) 
    error_estimate_sampled = data.frame(error_estimate_sampled)
    error_estimate_sampled$dividedVar = error_estimate$meanVar/3
    
    ## transcriptome distance measurment, within category, sampled ##
    TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
      filter( sp1 != sp2 ) %>%
      rowwise() %>%
      group_by_all() %>%
      summarise(var_divergence = 
                  var(sampled.data[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] - 
                        sampled.data[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) - 
                  error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"] -
                  error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"]) %>%
      ungroup %>%
      spread( sp2, var_divergence, fill=0 ) %>%
      column_to_rownames( "sp1" ) %>%
      as.matrix
    
    TR_distances_with_errors_normalized = TR_distances_with_errors
    for(sp1 in organisms_list) {
      for(sp2 in organisms_list) {
        denominator = variance_estimates[variance_estimates$species == sp1 & variance_estimates$tissue == tissue, "TR_var"] + 
          variance_estimates[variance_estimates$species == sp2 & variance_estimates$tissue == tissue, "TR_var"] -
          error_estimate[error_estimate$species == sp1 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "dividedVar"] -
          error_estimate[error_estimate$species == sp2 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "dividedVar"] 
        TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/(denominator/2)
      }
    }
    
    TR_tree = nj( TR_distances_with_errors_normalized )
    TR_lengths_bootstrapped[i,gene_category] = sum(TR_tree$edge.length)

    ## translatome distance measurment, within category, sampled ##
    
    RP_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
      filter( sp1 != sp2 ) %>%
      rowwise() %>%
      group_by_all() %>%
      summarise(var_divergence = 
                  var(sampled.data[ , paste("median_",sp1,"_",tissue,".","RP",sep = "")] - 
                        sampled.data[ , paste("median_",sp2,"_",tissue,".","RP",sep = "")], na.rm=TRUE ) - 
                  error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "RP", "dividedVar"] -
                  error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "RP", "dividedVar"]) %>%
      ungroup %>%
      spread( sp2, var_divergence, fill=0 ) %>%
      column_to_rownames( "sp1" ) %>%
      as.matrix
    
    RP_distances_with_errors_normalized = RP_distances_with_errors
    for(sp1 in organisms_list) {
      for(sp2 in organisms_list) {
        denominator = variance_estimates[variance_estimates$species == sp1 & variance_estimates$tissue == tissue, "RP_var"] + 
          variance_estimates[variance_estimates$species == sp2 & variance_estimates$tissue == tissue, "RP_var"] -
          error_estimate[error_estimate$species == sp1 & error_estimate$type == "RP" & error_estimate$tissue == tissue, "dividedVar"] -
          error_estimate[error_estimate$species == sp2 & error_estimate$type == "RP" & error_estimate$tissue == tissue, "dividedVar"] 
        RP_distances_with_errors_normalized[sp1,sp2] = RP_distances_with_errors[sp1,sp2]/(denominator/2)
      }
    }
    
    RP_tree = nj( RP_distances_with_errors_normalized )
    RP_lengths_bootstrapped[i,gene_category] = sum(RP_tree$edge.length)
    
  }
}

