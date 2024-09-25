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

#Monotherapy----------------------------------------------------------------------
v1v2_pheno2<-v1v2_pheno1 %>% filter(Treatment_Type == "Immunotherapy only")

##PFS----------------------------------------------------------------------------
###Continuous--------------------------------------------------------------------
Mono_OC<-matrix(NA,45,3)
Mono_OC[,1]<-colnames(v1v2_pheno2)[c(8:19,30:62)]
rownames(Mono_OC)<-Mono_OC[,1]
for (i in colnames(v1v2_pheno2)[c(8:19,30:62)]) {
  #m<-median(T0_pheno2[,i],na.rm = T)
  #T0_pheno2$group<-ifelse(T0_pheno2[,i] >m , "Yes", "No")
  
  Mono_OC[i,2:3]<-round(as.numeric(summary(coxph(Surv(PFSDays, PFSEvent == 1) ~ v1v2_pheno2[,i]   + Sex + Age + Smoking, data=v1v2_pheno2))$coefficients[1,c(2,5)]),4)
  
}

Mono_OC<-Mono_OC[,-1]
colnames(Mono_OC)<-c("HR","p-value")  

#cox.zph(coxph(Surv(PFSurvivalDays, PFSEvent == 1) ~ v1v2_pheno2$gMDSC_Score + Sex + Age + Smoking, data=v1v2_pheno1))
write.csv(Mono_OC, file = "Mono_PFS_continuous_output.csv")

###Binary-------------------------------------------------------------------------
zeros_percentage <- colMeans(v1v2_pheno2[,c(8:19,30:62)]== 0, na.rm = TRUE) * 100  # Calculate percentage of zeros in each column
columns_over_80_percent_zeros <- names(zeros_percentage[zeros_percentage > 80])

v1v2_pheno3<-v1v2_pheno2[, !(names(v1v2_pheno2) %in% columns_over_80_percent_zeros)]

Mono_PB<-matrix(NA,38,6)
Mono_PB[,1]<-colnames(v1v2_pheno3)[c(8:17,28:55)]
rownames(Mono_PB)<-Mono_PB[,1]
for (i in colnames(v1v2_pheno3)[c(8:17,28:55)]) {
  m<-surv_cutpoint(v1v2_pheno3, time = "PFSDays", event = "PFSEvent", i, minprop = 0.2)$cutpoint["cutpoint"]$cutpoint
  v1v2_pheno3$group<-ifelse(v1v2_pheno3[,i] >m , "Yes", "No")
  Mono_PB[i,6]<-round(m,2)
  Mono_PB[i,5]<-round(as.numeric(summary(coxph(Surv(PFSDays, PFSEvent == 1) ~ group   + Sex + Age + Smoking, data=v1v2_pheno3))$coefficients[1,5]),4)
  Mono_PB[i,2:4]<-round(as.numeric(summary(coxph(Surv(PFSDays, PFSEvent == 1) ~ group   + Sex + Age + Smoking, data=v1v2_pheno3))$conf.int[1,c(1,3,4)]),4)
}

Mono_PB<-Mono_PB[,-1]
colnames(Mono_PB)<-c("HR","Lower","Higher","p_value","cut-off")  


write.csv(Mono_PB, file = "Mono_PFS_binary_output.csv")

###Forest Plot-------------------------------------------------------------------
Mono_PB<-as.data.frame(Mono_PB)
Mono_PB<-Mono_PB[order(Mono_PB$p_value),]
Mono_PB<-Mono_PB[order(Mono_PB$HR, decreasing = T),]
PFS_Sig<-Mono_PB %>% filter(p_value < 0.05)

orText <- cbind(c("Immune Variable",rownames(PFS_Sig)),
                c('Hazard Ratio (95% CI)',paste(format(round(as.numeric(PFS_Sig$HR),2),nsmall = 2),' (',format(round(as.numeric(PFS_Sig$Lower),2),
                                                                                                               nsmall = 2),', ',
                                                format(round(as.numeric(PFS_Sig$Higher,2)),nsmall = 2),')',sep = '')),
                c("Optimal Cut-off",PFS_Sig$`cut-off`))


library(forestplot)
forestplot(labeltext = orText, 
           hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:4, col = "#000044"),
                             "12" = gpar(lty = 1,columns = 1:4, col = "#000044")),
           #graph.pos = 2, 
           
           mean = c(NA,as.numeric(PFS_Sig$HR)),
           lower = c(NA,as.numeric(PFS_Sig$Lower)),
           upper = c(NA,as.numeric(PFS_Sig$Higher)),
           #title="IPS DMCs NR3C1 Enrichment",
           #txt_gp=fpTxtGp(label=gpar(cex=1.3),ticks=gpar(cex=1.3)),
           is.summary=c(TRUE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE, FALSE, FALSE,FALSE,FALSE,FALSE,FALSE),
           txt_gp = fpTxtGp(summary = list(
             gpar(cex=1.2)),
             ticks=gpar(cex=1.2),xlab=gpar(cex=1.2)
           )
           , xlab = "Hazard Ratio",
           
           # xticks = c(0,1,2,3,4,5,6,7,8,9),
           xlog = T,
           col=fpColors(box = "royalblue",
                        line = "darkblue")
           ,
           fn.ci_norm = fpDrawPointCI, pch = 16, cex=2, 
           zero = 1
           , 
           lineheight = "auto"
           ,
           boxsize=0.15
           , 
           colgap=unit(8,"mm")
)


##OS-----------------------------------------------------------------------------
###Continuous--------------------------------------------------------------------
Mono_OC<-matrix(NA,45,3)
Mono_OC[,1]<-colnames(v1v2_pheno2)[c(8:19,30:62)]
rownames(Mono_OC)<-Mono_OC[,1]
for (i in colnames(v1v2_pheno2)[c(8:19,30:62)]) {
  #m<-median(T0_pheno2[,i],na.rm = T)
  #T0_pheno2$group<-ifelse(T0_pheno2[,i] >m , "Yes", "No")
  
  Mono_OC[i,2:3]<-round(as.numeric(summary(coxph(Surv(OSDays, OSEvent == 1) ~ v1v2_pheno2[,i]   + Sex + Age + Smoking, data=v1v2_pheno2))$coefficients[1,c(2,5)]),4)
  
}

Mono_OC<-Mono_OC[,-1]
colnames(Mono_OC)<-c("HR","p-value")  

#cox.zph(coxph(Surv(PFSurvivalDays, OSEvent == 1) ~ v1v2_pheno2$gMDSC_Score + Sex + Age + Smoking, data=v1v2_pheno1))
write.csv(Mono_OC, file = "Mono_OS_continuous_output.csv")

###Binary-------------------------------------------------------------------------
zeros_percentage <- colMeans(v1v2_pheno2[,c(8:19,30:62)]== 0, na.rm = TRUE) * 100  # Calculate percentage of zeros in each column
columns_over_80_percent_zeros <- names(zeros_percentage[zeros_percentage > 80])

v1v2_pheno3<-v1v2_pheno2[, !(names(v1v2_pheno2) %in% columns_over_80_percent_zeros)]

Mono_OB<-matrix(NA,38,6)
Mono_OB[,1]<-colnames(v1v2_pheno3)[c(8:17,28:55)]
rownames(Mono_OB)<-Mono_OB[,1]
for (i in colnames(v1v2_pheno3)[c(8:17,28:55)]) {
  m<-surv_cutpoint(v1v2_pheno3, time = "OSDays", event = "OSEvent", i, minprop = 0.2)$cutpoint["cutpoint"]$cutpoint
  v1v2_pheno3$group<-ifelse(v1v2_pheno3[,i] >m , "Yes", "No")
  Mono_OB[i,6]<-round(m,2)
  Mono_OB[i,5]<-round(as.numeric(summary(coxph(Surv(OSDays, OSEvent == 1) ~ group   + Sex + Age + Smoking, data=v1v2_pheno3))$coefficients[1,5]),4)
  Mono_OB[i,2:4]<-round(as.numeric(summary(coxph(Surv(OSDays, OSEvent == 1) ~ group   + Sex + Age + Smoking, data=v1v2_pheno3))$conf.int[1,c(1,3,4)]),4)
}

Mono_OB<-Mono_OB[,-1]
colnames(Mono_OB)<-c("HR","Lower","Higher","p_value","cut-off")  


write.csv(Mono_OB, file = "Mono_OS_binary_output.csv")

###Forest Plot-------------------------------------------------------------------
Mono_OB<-as.data.frame(Mono_OB)
Mono_OB<-Mono_OB[order(Mono_OB$p_value),]
Mono_OB<-Mono_OB[order(Mono_OB$HR, decreasing = T),]
OS_Sig<-Mono_OB %>% filter(p_value < 0.05)

orText <- cbind(c("Immune Variable",rownames(OS_Sig)),
                c('Hazard Ratio (95% CI)',paste(format(round(as.numeric(OS_Sig$HR),2),nsmall = 2),' (',format(round(as.numeric(OS_Sig$Lower),2),
                                                                                                               nsmall = 2),', ',
                                                format(round(as.numeric(OS_Sig$Higher,2)),nsmall = 2),')',sep = '')),
                c("Optimal Cut-off",OS_Sig$`cut-off`))



forestplot(labeltext = orText, 
           hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:4, col = "#000044"),
                             "14" = gpar(lty = 1,columns = 1:4, col = "#000044")),
           #graph.pos = 2, 
           
           mean = c(NA,as.numeric(OS_Sig$HR)),
           lower = c(NA,as.numeric(OS_Sig$Lower)),
           upper = c(NA,as.numeric(OS_Sig$Higher)),
           #title="IPS DMCs NR3C1 Enrichment",
           #txt_gp=fpTxtGp(label=gpar(cex=1.3),ticks=gpar(cex=1.3)),
           is.summary=c(TRUE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE, FALSE, FALSE,FALSE,FALSE,FALSE,FALSE),
           txt_gp = fpTxtGp(summary = list(
             gpar(cex=1.2)),
             ticks=gpar(cex=1.2),xlab=gpar(cex=1.2)
           )
           , xlab = "Hazard Ratio",
           
           # xticks = c(0,1,2,3,4,5,6,7,8,9),
           xlog = T,
           col=fpColors(box = "royalblue",
                        line = "darkblue")
           ,
           fn.ci_norm = fpDrawPointCI, pch = 16, cex=2, 
           zero = 1
           , 
           lineheight = "auto"
           ,
           boxsize=0.15
           , 
           colgap=unit(8,"mm")
)




#Multitherapy----------------------------------------------------------------------
v1v2_pheno2<-v1v2_pheno1 

##PFS----------------------------------------------------------------------------
###Continuous--------------------------------------------------------------------
Mono_OC<-matrix(NA,45,3)
Mono_OC[,1]<-colnames(v1v2_pheno2)[c(8:19,30:62)]
rownames(Mono_OC)<-Mono_OC[,1]
for (i in colnames(v1v2_pheno2)[c(8:19,30:62)]) {
  #m<-median(T0_pheno2[,i],na.rm = T)
  #T0_pheno2$group<-ifelse(T0_pheno2[,i] >m , "Yes", "No")
  
  Mono_OC[i,2:3]<-round(as.numeric(summary(coxph(Surv(PFSDays, PFSEvent == 1) ~ v1v2_pheno2[,i]   + Sex + Age + Smoking, data=v1v2_pheno2))$coefficients[1,c(2,5)]),4)
  
}

Mono_OC<-Mono_OC[,-1]
colnames(Mono_OC)<-c("HR","p-value")  

#cox.zph(coxph(Surv(PFSurvivalDays, PFSEvent == 1) ~ v1v2_pheno2$gMDSC_Score + Sex + Age + Smoking, data=v1v2_pheno1))
write.csv(Mono_OC, file = "Multi_PFS_continuous_output.csv")

###Binary-------------------------------------------------------------------------
zeros_percentage <- colMeans(v1v2_pheno2[,c(8:19,30:62)]== 0, na.rm = TRUE) * 100  # Calculate percentage of zeros in each column
columns_over_80_percent_zeros <- names(zeros_percentage[zeros_percentage > 80])

v1v2_pheno3<-v1v2_pheno2[, !(names(v1v2_pheno2) %in% columns_over_80_percent_zeros)]

Mono_PB<-matrix(NA,38,6)
Mono_PB[,1]<-colnames(v1v2_pheno3)[c(8:17,28:55)]
rownames(Mono_PB)<-Mono_PB[,1]
for (i in colnames(v1v2_pheno3)[c(8:17,28:55)]) {
  m<-surv_cutpoint(v1v2_pheno3, time = "PFSDays", event = "PFSEvent", i, minprop = 0.2)$cutpoint["cutpoint"]$cutpoint
  v1v2_pheno3$group<-ifelse(v1v2_pheno3[,i] >m , "Yes", "No")
  Mono_PB[i,6]<-round(m,2)
  Mono_PB[i,5]<-round(as.numeric(summary(coxph(Surv(PFSDays, PFSEvent == 1) ~ group   + Sex + Age + Smoking, data=v1v2_pheno3))$coefficients[1,5]),4)
  Mono_PB[i,2:4]<-round(as.numeric(summary(coxph(Surv(PFSDays, PFSEvent == 1) ~ group   + Sex + Age + Smoking, data=v1v2_pheno3))$conf.int[1,c(1,3,4)]),4)
}

Mono_PB<-Mono_PB[,-1]
colnames(Mono_PB)<-c("HR","Lower","Higher","p_value","cut-off")  


write.csv(Mono_PB, file = "Multi_PFS_binary_output.csv")

###Forest Plot-------------------------------------------------------------------
Mono_PB<-as.data.frame(Mono_PB)
Mono_PB<-Mono_PB[order(Mono_PB$p_value),]
Mono_PB<-Mono_PB[order(Mono_PB$HR, decreasing = T),]
PFS_Sig<-Mono_PB %>% filter(p_value < 0.05)

orText <- cbind(c("Immune Variable",rownames(PFS_Sig)),
                c('Hazard Ratio (95% CI)',paste(format(round(as.numeric(PFS_Sig$HR),2),nsmall = 2),' (',format(round(as.numeric(PFS_Sig$Lower),2),
                                                                                                               nsmall = 2),', ',
                                                format(round(as.numeric(PFS_Sig$Higher,2)),nsmall = 2),')',sep = '')),
                c("Optimal Cut-off",PFS_Sig$`cut-off`))


library(forestplot)
forestplot(labeltext = orText, 
           hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:4, col = "#000044"),
                             "9" = gpar(lty = 1,columns = 1:4, col = "#000044")),
           #graph.pos = 2, 
           
           mean = c(NA,as.numeric(PFS_Sig$HR)),
           lower = c(NA,as.numeric(PFS_Sig$Lower)),
           upper = c(NA,as.numeric(PFS_Sig$Higher)),
           #title="IPS DMCs NR3C1 Enrichment",
           #txt_gp=fpTxtGp(label=gpar(cex=1.3),ticks=gpar(cex=1.3)),
           is.summary=c(TRUE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE, FALSE, FALSE,FALSE,FALSE,FALSE,FALSE),
           txt_gp = fpTxtGp(summary = list(
             gpar(cex=1.2)),
             ticks=gpar(cex=1.2),xlab=gpar(cex=1.2)
           )
           , xlab = "Hazard Ratio",
           
           # xticks = c(0,1,2,3,4,5,6,7,8,9),
           xlog = T,
           col=fpColors(box = "royalblue",
                        line = "darkblue")
           ,
           fn.ci_norm = fpDrawPointCI, pch = 16, cex=2, 
           zero = 1
           , 
           lineheight = "auto"
           ,
           boxsize=0.15
           , 
           colgap=unit(8,"mm")
)


##OS-----------------------------------------------------------------------------
###Continuous--------------------------------------------------------------------
Mono_OC<-matrix(NA,45,3)
Mono_OC[,1]<-colnames(v1v2_pheno2)[c(8:19,30:62)]
rownames(Mono_OC)<-Mono_OC[,1]
for (i in colnames(v1v2_pheno2)[c(8:19,30:62)]) {
  #m<-median(T0_pheno2[,i],na.rm = T)
  #T0_pheno2$group<-ifelse(T0_pheno2[,i] >m , "Yes", "No")
  
  Mono_OC[i,2:3]<-round(as.numeric(summary(coxph(Surv(OSDays, OSEvent == 1) ~ v1v2_pheno2[,i]   + Sex + Age + Smoking, data=v1v2_pheno2))$coefficients[1,c(2,5)]),4)
  
}

Mono_OC<-Mono_OC[,-1]
colnames(Mono_OC)<-c("HR","p-value")  

#cox.zph(coxph(Surv(PFSurvivalDays, OSEvent == 1) ~ v1v2_pheno2$gMDSC_Score + Sex + Age + Smoking, data=v1v2_pheno1))
write.csv(Mono_OC, file = "Multi_OS_continuous_output.csv")

###Binary-------------------------------------------------------------------------
zeros_percentage <- colMeans(v1v2_pheno2[,c(8:19,30:62)]== 0, na.rm = TRUE) * 100  # Calculate percentage of zeros in each column
columns_over_80_percent_zeros <- names(zeros_percentage[zeros_percentage > 80])

v1v2_pheno3<-v1v2_pheno2[, !(names(v1v2_pheno2) %in% columns_over_80_percent_zeros)]

Mono_OB<-matrix(NA,38,6)
Mono_OB[,1]<-colnames(v1v2_pheno3)[c(8:17,28:55)]
rownames(Mono_OB)<-Mono_OB[,1]
for (i in colnames(v1v2_pheno3)[c(8:17,28:55)]) {
  m<-surv_cutpoint(v1v2_pheno3, time = "OSDays", event = "OSEvent", i, minprop = 0.2)$cutpoint["cutpoint"]$cutpoint
  v1v2_pheno3$group<-ifelse(v1v2_pheno3[,i] >m , "Yes", "No")
  Mono_OB[i,6]<-round(m,2)
  Mono_OB[i,5]<-round(as.numeric(summary(coxph(Surv(OSDays, OSEvent == 1) ~ group   + Sex + Age + Smoking, data=v1v2_pheno3))$coefficients[1,5]),4)
  Mono_OB[i,2:4]<-round(as.numeric(summary(coxph(Surv(OSDays, OSEvent == 1) ~ group   + Sex + Age + Smoking, data=v1v2_pheno3))$conf.int[1,c(1,3,4)]),4)
}

Mono_OB<-Mono_OB[,-1]
colnames(Mono_OB)<-c("HR","Lower","Higher","p_value","cut-off")  


write.csv(Mono_OB, file = "Multi_OS_binary_output.csv")

###Forest Plot-------------------------------------------------------------------
Mono_OB<-as.data.frame(Mono_OB)
Mono_OB<-Mono_OB[order(Mono_OB$p_value),]
Mono_OB<-Mono_OB[order(Mono_OB$HR, decreasing = T),]
OS_Sig<-Mono_OB %>% filter(p_value < 0.05)

temp_row <- OS_Sig[4, ]  # Store the values of the 4th row temporarily
OS_Sig[4, ] <- OS_Sig[5, ]   # Assign the values of the 5th row to the 4th row
OS_Sig[5, ] <- temp_row
rownames(OS_Sig)[4:5]<-c("Eos","Ext_LMLratio")

orText <- cbind(c("Immune Variable",rownames(OS_Sig)),
                c('Hazard Ratio (95% CI)',paste(format(round(as.numeric(OS_Sig$HR),2),nsmall = 2),' (',format(round(as.numeric(OS_Sig$Lower),2),
                                                                                                              nsmall = 2),', ',
                                                format(round(as.numeric(OS_Sig$Higher,2)),nsmall = 2),')',sep = '')),
                c("Optimal Cut-off",OS_Sig$`cut-off`))



forestplot(labeltext = orText, 
           hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:4, col = "#000044"),
                             "11" = gpar(lty = 1,columns = 1:4, col = "#000044")),
           #graph.pos = 2, 
           
           mean = c(NA,as.numeric(OS_Sig$HR)),
           lower = c(NA,as.numeric(OS_Sig$Lower)),
           upper = c(NA,as.numeric(OS_Sig$Higher)),
           #title="IPS DMCs NR3C1 Enrichment",
           #txt_gp=fpTxtGp(label=gpar(cex=1.3),ticks=gpar(cex=1.3)),
           is.summary=c(TRUE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE, FALSE, FALSE,FALSE,FALSE,FALSE,FALSE),
           txt_gp = fpTxtGp(summary = list(
             gpar(cex=1.2)),
             ticks=gpar(cex=1.2),xlab=gpar(cex=1.2)
           )
           , xlab = "Hazard Ratio",
           
           # xticks = c(0,1,2,3,4,5,6,7,8,9),
           xlog = T,
           col=fpColors(box = "royalblue",
                        line = "darkblue")
           ,
           fn.ci_norm = fpDrawPointCI, pch = 16, cex=2, 
           zero = 1
           , 
           lineheight = "auto"
           ,
           boxsize=0.15
           , 
           colgap=unit(8,"mm")
)


