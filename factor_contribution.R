library(relaimpo)
library(plyr)

myPallet = c("#3399CC","#339900","#FF6600")
names(myPallet) = c("Brain", "Liver", "Testis")

#### load expression data ####
expression_filtered = readRDS("Data/expression_vs_factors.rds")

#### Brain contributions ####
fit_TR_Brain = lm(TR_var_Brain ~ pLI + HI + Tau_Brain + TR_median_Brain + Age + Ohno + PrC + maSigPro_Brain + dnds, data = expression_filtered)
summary(fit_TR_Brain)


TR_Brain_contributions = calc.relimp(fit_TR_Brain,type=c("lmg","last","first","pratt"),
            rela=TRUE)@lmg

names(TR_Brain_contributions) = c("pLI","HI","Tau","expression","Age","Ohno","PrC","maSigPro","omega")
TR_Brain2merge = data.frame(matrix(c(names(TR_Brain_contributions),TR_Brain_contributions,rep("Brain",length(TR_Brain_contributions)),rep("TR",length(TR_Brain_contributions))),ncol=4))

fit_RP_Brain = lm(RP_var_Brain ~ pLI + HI + Tau_Brain + RP_median_Brain + Age + Ohno + PrC + maSigPro_Brain + dnds, data = expression_filtered)
summary(fit_RP_Brain)
RP_Brain_contributions = calc.relimp(fit_RP_Brain,type=c("lmg","last","first","pratt"),
            rela=TRUE)@lmg
names(RP_Brain_contributions) = c("pLI","HI","Tau","expression","Age","Ohno","PrC","maSigPro","omega")
RP_Brain2merge = data.frame(matrix(c(names(RP_Brain_contributions),RP_Brain_contributions,rep("Brain",length(RP_Brain_contributions)),rep("RP",length(RP_Brain_contributions))),ncol=4))

#### Liver contributions ####
fit_TR_Liver = lm(TR_var_Liver ~ pLI + HI + Tau_Liver + TR_median_Liver + Age + Ohno + PrC + maSigPro_Liver + dnds, data = expression_filtered)
summary(fit_TR_Liver)
TR_Liver_contributions = calc.relimp(fit_TR_Liver,type=c("lmg","last","first","pratt"),
            rela=TRUE)@lmg
names(TR_Liver_contributions) = c("pLI","HI","Tau","expression","Age","Ohno","PrC","maSigPro","omega")
TR_Liver2merge = data.frame(matrix(c(names(TR_Liver_contributions),TR_Liver_contributions,rep("Liver",length(TR_Liver_contributions)),rep("TR",length(TR_Liver_contributions))),ncol=4))

fit_RP_Liver = lm(RP_var_Liver ~ pLI + HI + Tau_Liver + RP_median_Liver + Age + Ohno + PrC + maSigPro_Liver + dnds, data = expression_filtered)
summary(fit_TR_Liver)
RP_Liver_contributions = calc.relimp(fit_RP_Liver,type=c("lmg","last","first","pratt"),
            rela=TRUE)@lmg
names(RP_Liver_contributions) = c("pLI","HI","Tau","expression","Age","Ohno","PrC","maSigPro","omega")
RP_Liver2merge = data.frame(matrix(c(names(RP_Liver_contributions),RP_Liver_contributions,rep("Liver",length(RP_Liver_contributions)),rep("RP",length(RP_Liver_contributions))),ncol=4))

#### Testis contributions ####
fit_TR_Testis = lm(TR_var_Testis ~ pLI + HI + Tau_Testis + TR_median_Testis + Age + Ohno + PrC + maSigPro_Testis + dnds, data = expression_filtered)
summary(fit_TR_Testis)
TR_Testis_contributions = calc.relimp(fit_TR_Testis,type=c("lmg","last","first","pratt"),
                                     rela=TRUE)@lmg
names(TR_Testis_contributions) = c("pLI","HI","Tau","expression","Age","Ohno","PrC","maSigPro","omega")
TR_Testis2merge = data.frame(matrix(c(names(TR_Testis_contributions),TR_Testis_contributions,rep("Testis",length(TR_Testis_contributions)),rep("TR",length(TR_Testis_contributions))),ncol=4))

fit_RP_Testis = lm(RP_var_Testis ~ pLI + HI + Tau_Testis + RP_median_Testis + Age + Ohno + PrC + maSigPro_Testis + dnds, data = expression_filtered)
summary(fit_RP_Testis)
RP_Testis_contributions = calc.relimp(fit_RP_Testis,type=c("lmg","last","first","pratt"),
                                     rela=TRUE)@lmg
names(RP_Testis_contributions) = c("pLI","HI","Tau","expression","Age","Ohno","PrC","maSigPro","omega")
RP_Testis2merge = data.frame(matrix(c(names(RP_Testis_contributions),RP_Testis_contributions,rep("Testis",length(RP_Testis_contributions)),rep("RP",length(RP_Testis_contributions))),ncol=4))


#### plot factor contributions ####

contributions = rbind(TR_Brain2merge,RP_Brain2merge,TR_Liver2merge,RP_Liver2merge,TR_Testis2merge,RP_Testis2merge)
names(contributions) = c("feature","value","tissue","layer")

contributions$layer = as.vector(contributions$layer)

contributions[contributions$layer == "TR","layer"] = "transcriptome"
contributions[contributions$layer == "RP","layer"] = "translatome"

contributions$value = as.numeric(as.vector(contributions$value))
contributions$feature = factor(contributions$feature, levels = c("expression","Tau","omega","pLI","HI","Age","PrC","Ohno", "maSigPro"))
contributions$feature = revalue(contributions$feature, c("expression" = "expr","omega" = "dn/ds","maSigPro" = "dev"))


ggplot(contributions[contributions$feature %in% c("expr","Tau","dn/ds","pLI","HI","Age"),],aes(x = feature,y = value, fill = tissue,  alpha = layer)) +
  geom_bar(stat = "identity",position = "dodge",color = "black") +
  facet_grid(tissue ~ .) +
  scale_alpha_manual(values = c(0.2,1)) +
  scale_fill_manual(values = myPallet) +
#  scale_color_manual(values = myPallet) +
  theme_bw() +
  theme(panel.border =  element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  axis.line = element_line(colour = "black"), strip.background = element_rect(fill="white",size=2)) +
  #  theme(axis.text.x=element_text(angle=90,hjust=1,size = 22)) +
  theme(axis.text.x=element_text(size = 22)) +
  theme(axis.title = element_text(face = "bold",size = 24)) +
  theme(legend.text = element_text(size = 24)) +
  theme(legend.title = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 24)) +
  theme(axis.line.x = element_line(color = "white")) +
  theme(axis.ticks.x =  element_blank()) +
  theme(strip.text.y = element_text(size=24)) +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  labs(x = paste(""), y = paste("relative contributions")) +
  ylim(0,0.8) +
  scale_x_discrete(position = "top") 
