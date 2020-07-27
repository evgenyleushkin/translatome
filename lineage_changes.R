library(DESeq2)
library(IHW)
library(stringr)

organisms_list = c("Human","Macaque","Mouse","Opossum","Platypus","Chicken")

#### load expression and DESeq2 data ####

expression_filtered = read.table("Data/expression_filtered.txt",header = T)
lfc_filtered = read.table("Data/lfc_filtered.txt", sep = "\t", header = T)
lfcSE_filtered = read.table("Data/lfcSE_filtered.txt", sep = "\t", header = T)

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

####  variance estimates ####

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

##### blowing coefficients #####

blowing_coefficient = data.frame(matrix(NA,ncol = 4))
names(blowing_coefficient) = c("organism1","organism2","tissue","coefficient")
blowing_coefficient = blowing_coefficient[-1,]

for(tissue in c("Brain","Liver","Testis")) {
  for(organism1 in organisms_list) {
    for(organism2 in organisms_list) {
      new_coefficient = data.frame(matrix(c(organism1,organism2,tissue,(variance_estimates[variance_estimates$species ==  organism1& variance_estimates$tissue == tissue,"RP_var"] + 
                                                                          variance_estimates[variance_estimates$species == organism2 & variance_estimates$tissue == tissue,"RP_var"] - 
                                                                          error_estimate[error_estimate$species == organism1 & error_estimate$tissue == tissue & error_estimate$type == "RP","dividedVar"] -
                                                                          error_estimate[error_estimate$species == organism2 & error_estimate$tissue == tissue & error_estimate$type == "RP","dividedVar"])/
                                              (variance_estimates[variance_estimates$species == organism1 & variance_estimates$tissue == tissue,"TR_var"] + 
                                                 variance_estimates[variance_estimates$species == organism2 & variance_estimates$tissue == tissue,"TR_var"] - 
                                                 error_estimate[error_estimate$species == organism1 & error_estimate$tissue == tissue & error_estimate$type == "TR","dividedVar"] -
                                                 error_estimate[error_estimate$species == organism2 & error_estimate$tissue == tissue & error_estimate$type == "TR","dividedVar"])),ncol=4))
      names(new_coefficient) = c("organism1","organism2","tissue","coefficient")
      blowing_coefficient = rbind(blowing_coefficient,new_coefficient)
    }
  }
}

blowing_coefficient$coefficient = as.numeric(as.vector(blowing_coefficient$coefficient))


#### calculate p_values ####

calc_p_value = function(organism1, organism2, tissue) {
  z_score =   (abs(lfc_filtered[, paste(organism1, organism2, tolower(tissue), "TR", sep =
                                          ".")]) * sqrt(blowing_coefficient[blowing_coefficient$organism1 == organism1 &
                                                                              blowing_coefficient$organism2 == organism2 &
                                                                              blowing_coefficient$tissue == tissue, "coefficient"]) - abs(lfc_filtered[, paste(organism1, organism2, tolower(tissue), "RP", sep =
                                                                                                                                                                 ".")])) /
    sqrt((lfcSE_filtered[, paste(organism1, organism2, tolower(tissue), "TR", sep =
                                   ".")] * sqrt(blowing_coefficient[blowing_coefficient$organism1 == organism1 &
                                                                      blowing_coefficient$organism2 == organism2 &
                                                                      blowing_coefficient$tissue == tissue, "coefficient"]))^2 + (lfcSE_filtered[, paste(organism1, organism2, tolower(tissue), "RP", sep =
                                                                                                                                                           ".")])^2)
  return(pnorm(z_score))
}


#### Lineage specific changes ####
## Brain ##
# directional #
p_values = calc_p_value("Human","Mouse","Brain")
names(p_values) = lfc_filtered$geneID
polarized_p_values =  p_values[calc_p_value("Human","Opossum","Brain") < 0.05]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Human_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

p_values = calc_p_value("Macaque","Mouse","Brain")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Macaque","Opossum","Brain") < 0.05]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Macaque_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

p_values = calc_p_value("Mouse","Macaque","Brain")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Mouse","Opossum","Brain") < 0.05]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Mouse_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))


# buffered #

p_values = 1 - calc_p_value("Human","Mouse","Brain")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Human","Opossum","Brain") > 0.95]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Human_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

p_values = 1 - calc_p_value("Macaque","Mouse","Brain")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Macaque","Opossum","Brain") > 0.95]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Macaque_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

p_values = 1 - calc_p_value("Mouse","Macaque","Brain")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Mouse","Opossum","Brain") > 0.95]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Mouse_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

## Testis ##
# directional #
p_values = calc_p_value("Human","Mouse","Testis")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Human","Opossum","Testis") < 0.05]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Human_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

p_values = calc_p_value("Macaque","Mouse","Testis")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Macaque","Opossum","Testis") < 0.05]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Macaque_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

p_values = calc_p_value("Mouse","Macaque","Testis")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Mouse","Opossum","Testis") < 0.05]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Mouse_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

# buffered #

p_values = 1 - calc_p_value("Human","Mouse","Testis")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Human","Opossum","Testis") > 0.95]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Human_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

p_values = 1 - calc_p_value("Macaque","Mouse","Testis")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Macaque","Opossum","Brain") > 0.95]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Macaque_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

p_values = 1 - calc_p_value("Mouse","Macaque","Testis")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Mouse","Opossum","Testis") > 0.95]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Mouse_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

## Liver ##
# directional #

p_values = calc_p_value("Macaque","Mouse","Liver")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Macaque","Opossum","Liver") < 0.05]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Macaque_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

p_values = calc_p_value("Mouse","Macaque","Liver")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Mouse","Opossum","Liver") < 0.05]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Mouse_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

# buffered #
p_values = 1 - calc_p_value("Macaque","Mouse","Liver")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Macaque","Opossum","Liver") > 0.95]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Macaque_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

p_values = 1 - calc_p_value("Mouse","Macaque","Liver")
names(p_values) = lfc_filtered$geneID
polarized_p_values = p_values[calc_p_value("Mouse","Opossum","Liver") > 0.95]
sum(p.adjust(polarized_p_values,method = "BH") < 0.05,na.rm=T)
Mouse_changes = names(na.omit(polarized_p_values[p.adjust(polarized_p_values,method = "BH") < 0.05]))

