library(ape)
library(cowplot)
library(tidyverse)

#### load expression data ####

expression_filtered = read.table("Data/expression_filtered.txt",header = T)

#### estimate of within-species vriance and measurement error ####

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

## set tissue: Berain, Liver or Testis

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

#### transcriptome layer divergence calculations, corrected ####

TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
  filter( sp1 != sp2 ) %>%
  rowwise() %>%
  group_by_all() %>%
  summarise(var_divergence = 
              var(expression_filtered[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] - 
                    expression_filtered[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) - 
              error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"] -
              error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"]) %>%
  ungroup %>%
  spread( sp2, var_divergence, fill=0 ) %>%
  column_to_rownames( "sp1" ) %>%
  as.matrix

TR_distances_with_errors_normalized = TR_distances_with_errors

variance_sum = 0
error_sum = 0
for(sp in organisms_list) {
  variance_sum = variance_sum + variance_estimates[variance_estimates$species == sp & variance_estimates$tissue == tissue, "TR_var"] 
  error_sum = error_sum + error_estimate[error_estimate$species == sp & error_estimate$type == "TR" & error_estimate$tissue == tissue, "dividedVar"]
}
denumerator = (variance_sum - error_sum)/length(organisms_list)

for(sp1 in organisms_list) {
  for(sp2 in organisms_list) {
    TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/denumerator
  }
}

TR_tree = nj( TR_distances_with_errors_normalized )
TR_div = sum(TR_tree$edge.length)

#### translatome layer divergence calculations, corrected ####

RP_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
  filter( sp1 != sp2 ) %>%
  rowwise() %>%
  group_by_all() %>%
  summarise(var_divergence = 
              var(expression_filtered[ , paste("median_",sp1,"_",tissue,".","RP",sep = "")] - 
                    expression_filtered[ , paste("median_",sp2,"_",tissue,".","RP",sep = "")], na.rm=TRUE ) - 
              error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "RP", "dividedVar"] -
              error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "RP", "dividedVar"]) %>%
  ungroup %>%
  spread( sp2, var_divergence, fill=0 ) %>%
  column_to_rownames( "sp1" ) %>%
  as.matrix


RP_distances_with_errors_normalized = RP_distances_with_errors

variance_sum = 0
error_sum = 0
for(sp in organisms_list) {
  variance_sum = variance_sum + variance_estimates[variance_estimates$species == sp & variance_estimates$tissue == tissue, "RP_var"] 
  error_sum = error_sum + error_estimate[error_estimate$species == sp & error_estimate$type == "RP" & error_estimate$tissue == tissue, "dividedVar"]
}
denumerator = (variance_sum - error_sum)/length(organisms_list)


for(sp1 in organisms_list) {
  for(sp2 in organisms_list) {
    RP_distances_with_errors_normalized[sp1,sp2] = RP_distances_with_errors[sp1,sp2]/denumerator
  }
}

RP_tree = nj( RP_distances_with_errors_normalized )
RP_div = sum(RP_tree$edge.length)

######## trees bootstrap ########

#### transcriptome layer ####

bootstrap_support1 = c()
bootstrap_support2 = c()
bootstrap_support3 = c()
bootstrap_support4 = c()

for(i in 1:100) {
  
  expression_bootstrapped = expression_filtered[sample(nrow(expression_filtered),replace=T),]
  
  #### estimate of within-species vriance and measurement error ####
  
  error_estimate = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( grepl("_R", colnames(expression_bootstrapped) )) %>%
    separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
    group_by( species, tissue, type ) %>%
    summarise( medianVar = median( genefilter::rowVars( expression_bootstrapped[,colname,drop=FALSE] ), na.rm=TRUE ) )
  error_estimate = data.frame(error_estimate)
  
  ####  variance estimates ####
  
  variance_estimates_RP_TR_TE = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( str_detect( colname, "median" ) ) %>%
    filter( colname != "Gene_ID" ) %>%
    separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
    #  filter( tissue == tissue) %>%
    spread( type, colname ) %>%
    group_by_all() %>%
    summarise( RP_var = var( expression_bootstrapped[ , RP ], na.rm=TRUE ) ) %>%
    ungroup %>%
    group_by_all() %>%
    summarise( TR_var = var( expression_bootstrapped[ , TR ], na.rm=TRUE ) ) %>%
    ungroup %>%
    group_by_all() %>%
    summarise( TE_var = var( expression_bootstrapped[ , RP ] - expression_bootstrapped[ , TR ], na.rm=TRUE ) ) %>%
    ungroup
  variance_estimates = data.frame(variance_estimates_RP_TR_TE)

  #### expression divergence calculations ####
  
  TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
    filter( sp1 != sp2 ) %>%
    rowwise() %>%
    group_by_all() %>%
    summarise(var_divergence =
                var(expression_bootstrapped[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] -
                      expression_bootstrapped[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) -
                error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "TR", "medianVar"] -
                error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "TR", "medianVar"]) %>%
    ungroup %>%
    spread( sp2, var_divergence, fill=0 ) %>%
    column_to_rownames( "sp1" ) %>%
    as.matrix
  
  TR_distances_with_errors_normalized = TR_distances_with_errors
  for(sp1 in organisms_list) {
    for(sp2 in organisms_list) {
      denominator = variance_estimates[variance_estimates$species == sp1 & variance_estimates$tissue == tissue, "TR_var"] +
        variance_estimates[variance_estimates$species == sp2 & variance_estimates$tissue == tissue, "TR_var"] -
        error_estimate[error_estimate$species == sp1 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "medianVar"] -
        error_estimate[error_estimate$species == sp2 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "medianVar"]
      TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/(denominator/2)
    }
  }
  
  TR_tree = nj( TR_distances_with_errors_normalized )
  
  bootstrap_support1 = c(bootstrap_support1,prop.clades(TR_tree,control_tree)[1])
  bootstrap_support2 = c(bootstrap_support2,prop.clades(TR_tree,control_tree)[2])
  bootstrap_support3 = c(bootstrap_support3,prop.clades(TR_tree,control_tree)[3])
  bootstrap_support4 = c(bootstrap_support4,prop.clades(TR_tree,control_tree)[4])
  
}

#### translatome layer ####

bootstrap_support1 = c()
bootstrap_support2 = c()
bootstrap_support3 = c()
bootstrap_support4 = c()

for(i in 1:100) {
  expression_bootstrapped = expression_filtered[sample(nrow(expression_filtered),replace=T),]
  
  #### estimate of within-species vriance and measurement error ####
  
  error_estimate = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( grepl("_R", colnames(expression_bootstrapped) )) %>%
    separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
    group_by( species, tissue, type ) %>%
    summarise( medianVar = median( genefilter::rowVars( expression_bootstrapped[,colname,drop=FALSE] ), na.rm=TRUE ) )
  error_estimate = data.frame(error_estimate)
  
  ####  variance estimates ####
  
  variance_estimates_RP_TR_TE = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( str_detect( colname, "median" ) ) %>%
    filter( colname != "Gene_ID" ) %>%
    separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
    #  filter( tissue == tissue) %>%
    spread( type, colname ) %>%
    group_by_all() %>%
    summarise( RP_var = var( expression_bootstrapped[ , RP ], na.rm=TRUE ) ) %>%
    ungroup %>%
    group_by_all() %>%
    summarise( TR_var = var( expression_bootstrapped[ , TR ], na.rm=TRUE ) ) %>%
    ungroup %>%
    group_by_all() %>%
    summarise( TE_var = var( expression_bootstrapped[ , RP ] - expression_bootstrapped[ , TR ], na.rm=TRUE ) ) %>%
    ungroup
  variance_estimates = data.frame(variance_estimates_RP_TR_TE)
 
  #### expression divergence calculations ####
  
  RP_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
    filter( sp1 != sp2 ) %>%
    rowwise() %>%
    group_by_all() %>%
    summarise(var_divergence = 
                var(expression_bootstrapped[ , paste("median_",sp1,"_",tissue,".","RP",sep = "")] - 
                      expression_bootstrapped[ , paste("median_",sp2,"_",tissue,".","RP",sep = "")], na.rm=TRUE ) - 
                error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "RP", "medianVar"] -
                error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "RP", "medianVar"]) %>%
    ungroup %>%
    spread( sp2, var_divergence, fill=0 ) %>%
    column_to_rownames( "sp1" ) %>%
    as.matrix
  
  RP_distances_with_errors_normalized = RP_distances_with_errors
  for(sp1 in organisms_list) {
    for(sp2 in organisms_list) {
      denominator = variance_estimates[variance_estimates$species == sp1 & variance_estimates$tissue == tissue, "RP_var"] + 
        variance_estimates[variance_estimates$species == sp2 & variance_estimates$tissue == tissue, "RP_var"] -
        error_estimate[error_estimate$species == sp1 & error_estimate$type == "RP" & error_estimate$tissue == tissue, "medianVar"] -
        error_estimate[error_estimate$species == sp2 & error_estimate$type == "RP" & error_estimate$tissue == tissue, "medianVar"] 
      RP_distances_with_errors_normalized[sp1,sp2] = RP_distances_with_errors[sp1,sp2]/(denominator/2)
    }
  }
  
  RP_tree = nj( RP_distances_with_errors_normalized )
  
  bootstrap_support1 = c(bootstrap_support1,prop.clades(RP_tree,control_tree)[1])
  bootstrap_support2 = c(bootstrap_support2,prop.clades(RP_tree,control_tree)[2])
  bootstrap_support3 = c(bootstrap_support3,prop.clades(RP_tree,control_tree)[3])
  bootstrap_support4 = c(bootstrap_support3,prop.clades(RP_tree,control_tree)[4])
  
}



