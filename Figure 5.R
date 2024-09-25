library(dplyr)
library(tibble)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(ggbeeswarm)
library(pROC)
library(survival) 
library(survminer)

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/HNSCR01/Analysis_032324/v1v2_pheno.RDATA")
for (col in names(v1v2_pheno)) {
  v1v2_pheno[[col]][is.nan(v1v2_pheno[[col]])] <- NA
  v1v2_pheno[[col]][is.infinite(v1v2_pheno[[col]])] <- NA
}

v1v2_pheno1<-v1v2_pheno %>% filter(is.na(Manual_Count) == T | Manual_Count != "Yes")
v1v2_pheno2<-v1v2_pheno1 %>% filter(Treatment_Type == "Immunotherapy only")

TMB_pheno<-read.csv("Updated_TMB.csv")
TMB_pheno<-TMB_pheno[complete.cases(TMB_pheno), ]
colnames(TMB_pheno)[1]<-"Subject_ID"
v1v2TMB_pheno<-inner_join(v1v2_pheno2,TMB_pheno,by = "Subject_ID")

#TMB and PFS Interaction------------------------------------------------------------
TMBPFS_output<-matrix(NA,45,5)
TMBPFS_output[,1]<-colnames(v1v2TMB_pheno)[c(8:19,30:62)]
rownames(TMBPFS_output)<-TMBPFS_output[,1]
for (i in colnames(v1v2TMB_pheno)[c(8:19,30:62)]) {
  TMBPFS_output[i,5]<-round(as.numeric(summary(coxph(Surv(PFSDays, PFSEvent == 1) ~ v1v2TMB_pheno[,"TMB"]*v1v2TMB_pheno[,i]   + Sex + Age + Smoking, data=v1v2TMB_pheno))$coefficient[6,5]),4)
  TMBPFS_output[i,2:4]<-round(as.numeric(summary(coxph(Surv(PFSDays, PFSEvent == 1) ~ v1v2TMB_pheno[,"TMB"]*v1v2TMB_pheno[,i]   + Sex + Age + Smoking, data=v1v2TMB_pheno))$conf.int[6,c(1,3,4)]),4)
}

TMBPFS_output<-TMBPFS_output[,-1]
colnames(TMBPFS_output)<-c("HR","Lower","Higher","p_value")  
write.csv(TMBPFS_output,file = "Interaction_TMB_PFS.csv")

#TMB and OS Interaction------------------------------------------------------------
TMBOS_output<-matrix(NA,45,5)
TMBOS_output[,1]<-colnames(v1v2TMB_pheno)[c(8:19,30:62)]
rownames(TMBOS_output)<-TMBOS_output[,1]
for (i in colnames(v1v2TMB_pheno)[c(8:19,30:62)]) {
  TMBOS_output[i,5]<-round(as.numeric(summary(coxph(Surv(OSDays, OSEvent == 1) ~ v1v2TMB_pheno[,"TMB"]*v1v2TMB_pheno[,i]   + Sex + Age + Smoking, data=v1v2TMB_pheno))$coefficient[6,5]),4)
  TMBOS_output[i,2:4]<-round(as.numeric(summary(coxph(Surv(OSDays, OSEvent == 1) ~ v1v2TMB_pheno[,"TMB"]*v1v2TMB_pheno[,i]   + Sex + Age + Smoking, data=v1v2TMB_pheno))$conf.int[6,c(1,3,4)]),4)
}

TMBOS_output<-TMBOS_output[,-1]
colnames(TMBOS_output)<-c("HR","Lower","Higher","p_value")  
write.csv(TMBOS_output,file = "Interaction_TMB_OS.csv")

#Monocyte Proportion Stratification-------------------------------------------------
library(gtools)
v1v2TMB_pheno$Mono_Quant<-quantcut(v1v2TMB_pheno$Mono, q = 2)
table(v1v2TMB_pheno$Mono_Quan)
v1v2TMB_pheno_1<-v1v2TMB_pheno %>% filter(Mono_Quant == "[3.37,8.29]")
cox.zph(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_1$TMB >3) +Sex + Age+Smoking, data=v1v2TMB_pheno_1))
summary(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_1$TMB >3) +Sex + Age+Smoking, data=v1v2TMB_pheno_1))

model <- survfit(Surv(OSDays, OSEvent == 1) ~ TMB>3,
                 data = v1v2TMB_pheno_1)
summary(model)
p<-ggsurvplot(model,  
              pval = T, #conf.int = TRUE,
              risk.table = T,
              ggtheme = theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               legend.text = element_text(size=15, 
                                                          face="bold"),
                               plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
                               axis.text.x = element_text(size = 13, face = "bold" ),
                               axis.text.y = element_text(size = 13, face = "bold" ),
                               axis.title =  element_text( size = 20 , face="bold")),
              #legend.labs =  c("TMB < 3 (n=11)", "TMB > 3 (n=12)"), 
              #linetype = "twodash",
              break.time.by = 100,
              title = "TMB and OS in Monocyte < 8.02 (%)",
              ylab = "Survival Proportion",
              xlab = "Time (Days)")

p$plot+scale_colour_manual(values = c("#4472C4","#EC491C"))



v1v2TMB_pheno_2<-v1v2TMB_pheno %>% filter(Mono_Quant == "(8.29,16.3]")
cox.zph(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_2$TMB > 3) +Sex + Age +Smoking, data=v1v2TMB_pheno_2))
summary(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_2$TMB > 3) +Sex + Age +Smoking, data=v1v2TMB_pheno_2))
model <- survfit(Surv(OSDays, OSEvent == 1) ~ TMB>3,
                 data = v1v2TMB_pheno_2)
summary(model)
p<-ggsurvplot(model,  
              pval = T, #conf.int = TRUE,
              risk.table = T,
              ggtheme = theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               legend.text = element_text(size=15, 
                                                          face="bold"),
                               plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
                               axis.text.x = element_text(size = 13, face = "bold" ),
                               axis.text.y = element_text(size = 13, face = "bold" ),
                               axis.title =  element_text( size = 20 , face="bold")),
              #legend.labs =  c("TMB < 3 (n=6)", "TMB > 3 (n=17)"), 
              #linetype = "twodash",
              break.time.by = 100,
              title = "HNSC TMB and OS in Monocyte > 8.02%",
              ylab = "Survival Proportion",
              xlab = "Time (Days)")

p$plot+scale_colour_manual(values = c("#4472C4","#EC491C"))







#CD4 Proportion Stratification-------------------------------------------------
v1v2TMB_pheno$Ext_CD4tot_Quant<-quantcut(v1v2TMB_pheno$Ext_CD4tot, q = 2)
table(v1v2TMB_pheno$Ext_CD4tot_Quan)
v1v2TMB_pheno_1<-v1v2TMB_pheno %>% filter(Ext_CD4tot_Quant == "[0.52,5.23]")
cox.zph(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_1$TMB >3) +Sex + Age+Smoking, data=v1v2TMB_pheno_1))
summary(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_1$TMB >3) +Sex + Age+Smoking, data=v1v2TMB_pheno_1))

model <- survfit(Surv(OSDays, OSEvent == 1) ~ TMB>3,
                 data = v1v2TMB_pheno_1)
summary(model)
p<-ggsurvplot(model,  
              pval = T, #conf.int = TRUE,
              risk.table = T,
              ggtheme = theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               legend.text = element_text(size=15, 
                                                          face="bold"),
                               plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
                               axis.text.x = element_text(size = 13, face = "bold" ),
                               axis.text.y = element_text(size = 13, face = "bold" ),
                               axis.title =  element_text( size = 20 , face="bold")),
              #legend.labs =  c("TMB < 3 (n=11)", "TMB > 3 (n=12)"), 
              #linetype = "twodash",
              break.time.by = 100,
              title = "TMB and OS in Ext_CD4totcyte < 8.02 (%)",
              ylab = "Survival Proportion",
              xlab = "Time (Days)")

p$plot+scale_colour_manual(values = c("#4472C4","#EC491C"))



v1v2TMB_pheno_2<-v1v2TMB_pheno %>% filter(Ext_CD4tot_Quant == "(5.23,17.8]")
cox.zph(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_2$TMB > 3) +Sex + Age +Smoking, data=v1v2TMB_pheno_2))
summary(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_2$TMB > 3) +Sex + Age +Smoking, data=v1v2TMB_pheno_2))
model <- survfit(Surv(OSDays, OSEvent == 1) ~ TMB>3,
                 data = v1v2TMB_pheno_2)
summary(model)
p<-ggsurvplot(model,  
              pval = T, #conf.int = TRUE,
              risk.table = T,
              ggtheme = theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               legend.text = element_text(size=15, 
                                                          face="bold"),
                               plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
                               axis.text.x = element_text(size = 13, face = "bold" ),
                               axis.text.y = element_text(size = 13, face = "bold" ),
                               axis.title =  element_text( size = 20 , face="bold")),
              #legend.labs =  c("TMB < 3 (n=6)", "TMB > 3 (n=17)"), 
              #linetype = "twodash",
              break.time.by = 100,
              title = "HNSC TMB and OS in Ext_CD4totcyte > 8.02%",
              ylab = "Survival Proportion",
              xlab = "Time (Days)")

p$plot+scale_colour_manual(values = c("#4472C4","#EC491C"))







#NK Count Stratification-------------------------------------------------
v1v2TMB_pheno$ExtAbsolute_NK_Quant<-quantcut(v1v2TMB_pheno$ExtAbsolute_NK, q = 2)
table(v1v2TMB_pheno$ExtAbsolute_NK_Quan)
v1v2TMB_pheno_1<-v1v2TMB_pheno %>% filter(ExtAbsolute_NK_Quant == "[0,134]")
cox.zph(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_1$TMB >3) +Sex + Age+Smoking, data=v1v2TMB_pheno_1))
summary(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_1$TMB >3) +Sex + Age+Smoking, data=v1v2TMB_pheno_1))

model <- survfit(Surv(OSDays, OSEvent == 1) ~ TMB>3,
                 data = v1v2TMB_pheno_1)
summary(model)
p<-ggsurvplot(model,  
              pval = T, #conf.int = TRUE,
              risk.table = T,
              ggtheme = theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               legend.text = element_text(size=15, 
                                                          face="bold"),
                               plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
                               axis.text.x = element_text(size = 13, face = "bold" ),
                               axis.text.y = element_text(size = 13, face = "bold" ),
                               axis.title =  element_text( size = 20 , face="bold")),
              #legend.labs =  c("TMB < 3 (n=11)", "TMB > 3 (n=12)"), 
              #linetype = "twodash",
              break.time.by = 100,
              title = "TMB and OS in ExtAbsolute_NKcyte < 8.02 (%)",
              ylab = "Survival Proportion",
              xlab = "Time (Days)")

p$plot+scale_colour_manual(values = c("#4472C4","#EC491C"))



v1v2TMB_pheno_2<-v1v2TMB_pheno %>% filter(ExtAbsolute_NK_Quant == "(134,448]")
cox.zph(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_2$TMB > 3) +Sex + Age +Smoking, data=v1v2TMB_pheno_2))
summary(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno_2$TMB > 3) +Sex + Age +Smoking, data=v1v2TMB_pheno_2))
model <- survfit(Surv(OSDays, OSEvent == 1) ~ TMB>3,
                 data = v1v2TMB_pheno_2)
summary(model)
p<-ggsurvplot(model,  
              pval = T, #conf.int = TRUE,
              risk.table = T,
              ggtheme = theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               legend.text = element_text(size=15, 
                                                          face="bold"),
                               plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
                               axis.text.x = element_text(size = 13, face = "bold" ),
                               axis.text.y = element_text(size = 13, face = "bold" ),
                               axis.title =  element_text( size = 20 , face="bold")),
              #legend.labs =  c("TMB < 3 (n=6)", "TMB > 3 (n=17)"), 
              #linetype = "twodash",
              break.time.by = 100,
              title = "HNSC TMB and OS in ExtAbsolute_NKcyte > 8.02%",
              ylab = "Survival Proportion",
              xlab = "Time (Days)")

p$plot+scale_colour_manual(values = c("#4472C4","#EC491C"))





