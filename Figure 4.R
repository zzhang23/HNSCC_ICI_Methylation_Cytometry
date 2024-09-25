library(dplyr)
library(tibble)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(ggbeeswarm)
library(pROC)
library(survival) 
library(survminer)

#TMB------------------------------------------------------------------------------
TMB_pheno<-read.csv("Updated_TMB.csv")
TMB_pheno<-TMB_pheno[complete.cases(TMB_pheno), ]

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/HNSCR01/Analysis_032324/v1v2_pheno.RDATA")
colnames(TMB_pheno)[1]<-"Subject_ID"
v1v2TMB_pheno<-inner_join(v1v2_pheno,TMB_pheno,by = "Subject_ID")

##Monotherapy------------------------------------------------------------------------
v1v2TMB_pheno1<-v1v2TMB_pheno %>% filter(Treatment_Type == "Immunotherapy only")

###PFS--------------------------------------------------------------------------------
#m<-surv_cutpoint(v1v2TMB_pheno1, time = "PFSDays", event = "PFSEvent", "TMB", minprop = 0.2)$cutpoint["cutpoint"]$cutpoint
m<-3
model <- survfit(Surv(PFSDays, PFSEvent == 1) ~ TMB>m,
                 data = v1v2TMB_pheno1)
summary(model)
p<-ggsurvplot(model,  
              # pval = T, #conf.int = TRUE,
              risk.table = T,
              ggtheme = theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               legend.text = element_text(size=15, 
                                                          face="bold"),
                               plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
                               axis.text.x = element_text(size = 13, face = "bold" ),
                               axis.text.y = element_text(size = 13, face = "bold" ),
                               axis.title =  element_text( size = 20 , face="bold")),
              legend.labs =  c("TMB <= 3", "TMB > 3"), 
              
              break.time.by = 100,
              title = "TMB and Progression-free Survival",
              ylab = "Survival Proportion",
              xlab = "Time (Days)")

p$plot+scale_colour_manual(values = c("#4472C4","#EC491C"))

cox.zph(coxph(Surv(PFSDays, PFSEvent == 1) ~ (v1v2TMB_pheno1$TMB> m) + Sex + Age + Smoking, data=v1v2TMB_pheno1))
summary(coxph(Surv(PFSDays, PFSEvent == 1) ~ (v1v2TMB_pheno1$TMB> m) + Sex + Age + Smoking, data=v1v2TMB_pheno1))

###OS----------------------------------------------------------------------------------
#m<-surv_cutpoint(v1v2TMB_pheno1, time = "OSDays", event = "OSEvent", "TMB", minprop = 0.2)$cutpoint["cutpoint"]$cutpoint
m<-3
model <- survfit(Surv(OSDays, OSEvent == 1) ~ TMB>m,
                 data = v1v2TMB_pheno1)
summary(model)
p<-ggsurvplot(model,  
              # pval = T, #conf.int = TRUE,
              risk.table = T,
              ggtheme = theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               legend.text = element_text(size=15, 
                                                          face="bold"),
                               plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
                               axis.text.x = element_text(size = 13, face = "bold" ),
                               axis.text.y = element_text(size = 13, face = "bold" ),
                               axis.title =  element_text( size = 20 , face="bold")),
              legend.labs =  c("TMB <= 3", "TMB > 3"), 
              
              break.time.by = 100,
              title = "HNSC TMB and OS",
              ylab = "Survival Proportion",
              xlab = "Time (Days)")

p$plot+scale_colour_manual(values = c("#4472C4","#EC491C"))


cox.zph(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno1$TMB> m) + Sex + Age + Smoking, data=v1v2TMB_pheno1))
summary(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2TMB_pheno1$TMB> m) + Sex + Age + Smoking, data=v1v2TMB_pheno1))

#PDL1 CPS--------------------------------------------------------------------------
CPS_pheno<-read.csv("Updated_PDL1_CPS.csv")
CPS_pheno<-CPS_pheno[complete.cases(CPS_pheno), ]

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/HNSCR01/Analysis_032324/v1v2_pheno.RDATA")
colnames(CPS_pheno)<-c("Subject_ID","PDL1_CPS")
v1v2PDL1_CPS_pheno<-inner_join(v1v2_pheno,CPS_pheno,by = "Subject_ID")


##Monotherapy------------------------------------------------------------------------
v1v2PDL1_CPS_pheno1<-v1v2PDL1_CPS_pheno %>% filter(Treatment_Type == "Immunotherapy only")

###PFS--------------------------------------------------------------------------------
m<-surv_cutpoint(v1v2PDL1_CPS_pheno1, time = "PFSDays", event = "PFSEvent", "PDL1_CPS", minprop = 0.2)$cutpoint["cutpoint"]$cutpoint
#m<-3
model <- survfit(Surv(PFSDays, PFSEvent == 1) ~ PDL1_CPS>m,
                 data = v1v2PDL1_CPS_pheno1)
summary(model)
p<-ggsurvplot(model,  
              # pval = T, #conf.int = TRUE,
              risk.table = T,
              ggtheme = theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               legend.text = element_text(size=15, 
                                                          face="bold"),
                               plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
                               axis.text.x = element_text(size = 13, face = "bold" ),
                               axis.text.y = element_text(size = 13, face = "bold" ),
                               axis.title =  element_text( size = 20 , face="bold")),
              legend.labs =  c("PDL1_CPS <= 3", "PDL1_CPS > 3"), 
              
              break.time.by = 100,
              title = "PDL1_CPS and Progression-free Survival",
              ylab = "Survival Proportion",
              xlab = "Time (Days)")

p$plot+scale_colour_manual(values = c("#4472C4","#EC491C"))

cox.zph(coxph(Surv(PFSDays, PFSEvent == 1) ~ (v1v2PDL1_CPS_pheno1$PDL1_CPS> m) + Sex + Age + Smoking, data=v1v2PDL1_CPS_pheno1))
summary(coxph(Surv(PFSDays, PFSEvent == 1) ~ (v1v2PDL1_CPS_pheno1$PDL1_CPS> m) + Sex + Age + Smoking, data=v1v2PDL1_CPS_pheno1))

###OS----------------------------------------------------------------------------------
m<-surv_cutpoint(v1v2PDL1_CPS_pheno1, time = "OSDays", event = "OSEvent", "PDL1_CPS", minprop = 0.2)$cutpoint["cutpoint"]$cutpoint
#m<-3
model <- survfit(Surv(OSDays, OSEvent == 1) ~ PDL1_CPS>m,
                 data = v1v2PDL1_CPS_pheno1)
summary(model)
p<-ggsurvplot(model,  
              # pval = T, #conf.int = TRUE,
              risk.table = T,
              ggtheme = theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               legend.text = element_text(size=15, 
                                                          face="bold"),
                               plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
                               axis.text.x = element_text(size = 13, face = "bold" ),
                               axis.text.y = element_text(size = 13, face = "bold" ),
                               axis.title =  element_text( size = 20 , face="bold")),
              legend.labs =  c("PDL1_CPS <= 3", "PDL1_CPS > 3"), 
              
              break.time.by = 100,
              title = "HNSC PDL1_CPS and OS",
              ylab = "Survival Proportion",
              xlab = "Time (Days)")

p$plot+scale_colour_manual(values = c("#4472C4","#EC491C"))


cox.zph(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2PDL1_CPS_pheno1$PDL1_CPS> m) + Sex + Age + Smoking, data=v1v2PDL1_CPS_pheno1))
summary(coxph(Surv(OSDays, OSEvent == 1) ~ (v1v2PDL1_CPS_pheno1$PDL1_CPS> m) + Sex + Age + Smoking, data=v1v2PDL1_CPS_pheno1))
