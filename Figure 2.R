library(dplyr)
library(tibble)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(ggbeeswarm)

load("v1v2_pheno.RDATA")
for (col in names(v1v2_pheno)) {
  v1v2_pheno[[col]][is.nan(v1v2_pheno[[col]])] <- NA
  v1v2_pheno[[col]][is.infinite(v1v2_pheno[[col]])] <- NA
}
rmsb<-v1v2_pheno %>% filter(PFSDays < 365 & PFSEvent ==0)
v1v2_pheno<-v1v2_pheno[!v1v2_pheno$Subject_ID %in% rmsb$Subject,]

v1v2_pheno$BenefitGroup<-ifelse(v1v2_pheno$PFSDays < 100, "Non-benefit",
                                ifelse(v1v2_pheno$PFSDays >= 365, "Durable benefit","Non-durable benefit"))


v1v2_pheno1<-v1v2_pheno %>% filter(is.na(Manual_Count) == T | Manual_Count != "Yes")

v1v2_pheno1<-v1v2_pheno1 %>% filter(Treatment_Type == "Immunotherapy only")

table(v1v2_pheno1$BenefitGroup)

v1v2_pheno1$BenefitGroup<-factor(v1v2_pheno1$BenefitGroup, levels=c("Durable benefit","Non-durable benefit","Non-benefit"))
v1v2_pheno2<-v1v2_pheno1%>% filter(BenefitGroup != "Non-benefit")
v1v2_pheno3<-v1v2_pheno1%>% filter(BenefitGroup != "Non-durable benefit")
v1v2_pheno4<-v1v2_pheno1%>% filter(BenefitGroup != "Durable benefit")

DBvsNDB_output<-matrix(NA,45,3)
DBvsNDB_output[,1]<-colnames(v1v2_pheno2)[c(8:19,30:62)]
rownames(DBvsNDB_output)<-DBvsNDB_output[,1]
for (i in colnames(v1v2_pheno2)[c(8:19,30:62)]) {
  DBvsNDB_output[i,2:3]<-round(as.numeric(summary(lm(as.formula(paste(i, "~", "BenefitGroup+Age+Sex+Smoking")), data = v1v2_pheno2))$coefficients[2,c(1,4)]),5)
}

DBvsNDB_output<-DBvsNDB_output[,-1]
colnames(DBvsNDB_output)<-c("coefficient","p-value")

DBvsNB_output<-matrix(NA,45,3)
DBvsNB_output[,1]<-colnames(v1v2_pheno3)[c(8:19,30:62)]
rownames(DBvsNB_output)<-DBvsNB_output[,1]
for (i in colnames(v1v2_pheno3)[c(8:19,30:62)]) {
  DBvsNB_output[i,2:3]<-round(as.numeric(summary(lm(as.formula(paste(i, "~", "BenefitGroup+Age+Sex+Smoking")), data = v1v2_pheno3))$coefficients[2,c(1,4)]),5)
}

DBvsNB_output<-DBvsNB_output[,-1]
colnames(DBvsNB_output)<-c("coefficient","p-value")

NBvsNDB_output<-matrix(NA,45,3)
NBvsNDB_output[,1]<-colnames(v1v2_pheno4)[c(8:19,30:62)]
rownames(NBvsNDB_output)<-NBvsNDB_output[,1]
for (i in colnames(v1v2_pheno4)[c(8:19,30:62)]) {
  NBvsNDB_output[i,2:3]<-round(as.numeric(summary(lm(as.formula(paste(i, "~", "BenefitGroup+Age+Sex+Smoking")), data = v1v2_pheno4))$coefficients[2,c(1,4)]),5)
}

NBvsNDB_output<-NBvsNDB_output[,-1]
colnames(NBvsNDB_output)<-c("coefficient","p-value")

#Making plots---------------------------------------------------------------------
# Define the comparisons for stat_compare_means
my_comparisons <- list(c("Durable benefit", "Non-benefit"), 
                       c("Durable benefit", "Non-durable benefit"), 
                       c("Non-durable benefit", "Non-benefit"))

# First plot
ggplot(v1v2_pheno1, aes(x=BenefitGroup, y=Ext_CD4nv_CD4mem, fill=BenefitGroup, label=Subject_ID)) + 
  geom_boxplot(fatten=1, outlier.shape=NA) + 
  geom_point(size=1.5) +
  stat_compare_means(label="p.format", method="t.test", size=4, comparisons=my_comparisons) +
  scale_fill_manual(values=c("#99CC66", "#FFCC66", "#E57373")) +  # Lighter colors
  ylab("Ratio CD4nv to CD4mem") +
  xlab("") +
  theme(legend.position="none", 
        legend.title=element_text(size=15, face="bold"), 
        legend.text=element_text(size=15, face="bold"), 
        axis.text.x=element_text(size=10, face="bold"), 
        axis.text.y=element_text(size=15, face="bold"), 
        axis.title=element_text(size=15, face="bold"), 
        plot.title=element_text(size=20, hjust=0.5)) + 
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(), 
        axis.line=element_line(colour="black"))

# Second plot
ggplot(v1v2_pheno1, aes(x=BenefitGroup, y=CD4nv, fill=BenefitGroup, label=Subject_ID)) + 
  geom_boxplot(fatten=1, outlier.shape=NA) + 
  geom_point(size=1.5) +
  stat_compare_means(label="p.format", method="t.test", size=4, comparisons=my_comparisons) +
  scale_fill_manual(values=c("#99CC66", "#FFCC66", "#E57373")) +  # Lighter colors
  ylab("Proportion CD4nv (%)") +
  xlab("") +
  theme(legend.position="none", 
        legend.title=element_text(size=15, face="bold"), 
        legend.text=element_text(size=15, face="bold"), 
        axis.text.x=element_text(size=10, face="bold"), 
        axis.text.y=element_text(size=15, face="bold"), 
        axis.title=element_text(size=15, face="bold"), 
        plot.title=element_text(size=20, hjust=0.5)) + 
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(), 
        axis.line=element_line(colour="black"))

# Third plot
ggplot(v1v2_pheno1, aes(x=BenefitGroup, y=Ext_CD4nv_pct, fill=BenefitGroup, label=Subject_ID)) + 
  geom_boxplot(fatten=1, outlier.shape=NA) + 
  geom_point(size=1.5) +
  stat_compare_means(label="p.format", method="t.test", size=4, comparisons=my_comparisons) +
  scale_fill_manual(values=c("#99CC66", "#FFCC66", "#E57373")) +  # Lighter colors
  ylab("Proportion CD4nv/CD4 (%)") +
  xlab("") +
  theme(legend.position="none", 
        legend.title=element_text(size=15, face="bold"), 
        legend.text=element_text(size=15, face="bold"), 
        axis.text.x=element_text(size=10, face="bold"), 
        axis.text.y=element_text(size=15, face="bold"), 
        axis.title=element_text(size=15, face="bold"), 
        plot.title=element_text(size=20, hjust=0.5)) + 
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(), 
        axis.line=element_line(colour="black"))
#Figure 2B-------------------------------------------------------------------------------------------

TMB_pheno<-read.csv("Updated_TMB.csv")
TMB_pheno<-TMB_pheno[complete.cases(TMB_pheno), ]

load("v1v2_pheno.RDATA")
colnames(TMB_pheno)[1]<-"Subject_ID"
v1v2TMB_pheno<-inner_join(v1v2_pheno,TMB_pheno,by = "Subject_ID")

rmsb<-v1v2TMB_pheno %>% filter(PFSDays < 365 & PFSEvent ==0)
v1v2TMB_pheno<-v1v2TMB_pheno[!v1v2TMB_pheno$Subject_ID %in% rmsb$Subject,]


v1v2TMB_pheno$BenefitGroup<-ifelse(v1v2TMB_pheno$PFSDays < 100, "Non-benefit",
                                   ifelse(v1v2TMB_pheno$PFSDays >= 365, "Durable benefit","Non-durable benefit"))

v1v2TMB_pheno1<-v1v2TMB_pheno %>% filter(Treatment_Type == "Immunotherapy only")

v1v2TMB_pheno1$BenefitGroup<-factor(v1v2TMB_pheno1$BenefitGroup, levels=c("Durable benefit","Non-durable benefit","Non-benefit"))

table(v1v2TMB_pheno1$BenefitGroup)

my_comparisons <- list( c("Durable benefit", "Non-benefit"), c("Durable benefit", "Non-durable benefit"), c("Non-durable benefit", "Non-benefit") )
ggplot(v1v2TMB_pheno1, aes(x=BenefitGroup, y=TMB,fill =BenefitGroup,label=Subject_ID)) + 
  geom_boxplot( fatten = 1,outlier.shape = NA)+geom_point(size = 1.5)+#geom_text(position =position_jitterdodge(seed = 1),size = 3)+
  stat_compare_means(label = "p.format",method = "t.test",size = 4,comparisons = my_comparisons)+#ggtitle("Immune Cell Projection in Multiple Sclerosis ")+
  scale_fill_manual(values=c("#99CC66", "#FFCC66", "#E57373"))+
  ylab("TMB /Mb")+xlab("")+
  theme( legend.position="none",legend.title = element_text( size=15, 
                                                             face="bold"),
         legend.text = element_text(size=15, 
                                    face="bold"),
         axis.text.x = element_text(size = 10, face = "bold" ),
         axis.text.y = element_text(size = 15, face = "bold" ),
         axis.title =  element_text( size = 15 , face="bold"),
         plot.title = element_text(size=20,hjust = 0.5))  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Figure 2C---------------------------------------------------------------------------

CPS_pheno<-read.csv("Updated_PDL1_CPS.csv")
CPS_pheno<-CPS_pheno[complete.cases(CPS_pheno), ]

load("v1v2_pheno.RDATA")
colnames(CPS_pheno)<-c("Subject_ID","PDL1_CPS")
v1v2CPS_pheno<-inner_join(v1v2_pheno,CPS_pheno,by = "Subject_ID")

rmsb<-v1v2CPS_pheno %>% filter(PFSDays < 365 & PFSEvent ==0)
v1v2CPS_pheno<-v1v2CPS_pheno[!v1v2CPS_pheno$Subject_ID %in% rmsb$Subject,]


v1v2CPS_pheno$BenefitGroup<-ifelse(v1v2CPS_pheno$PFSDays < 100, "Non-benefit",
                                   ifelse(v1v2CPS_pheno$PFSDays >= 365, "Durable benefit","Non-durable benefit"))

v1v2CPS_pheno1<-v1v2CPS_pheno %>% filter(Treatment_Type == "Immunotherapy only")

v1v2CPS_pheno1$BenefitGroup<-factor(v1v2CPS_pheno1$BenefitGroup, levels=c("Durable benefit","Non-durable benefit","Non-benefit"))
table(v1v2CPS_pheno1$BenefitGroup)
my_comparisons <- list( c("Durable benefit", "Non-benefit"), c("Durable benefit", "Non-durable benefit"), c("Non-durable benefit", "Non-benefit") )
ggplot(v1v2CPS_pheno1, aes(x=BenefitGroup, y=PDL1_CPS,fill =BenefitGroup,label=Subject_ID)) + 
  geom_boxplot( fatten = 1,outlier.shape = NA)+geom_point(size = 1.5)+#geom_text(position =position_jitterdodge(seed = 1),size = 3)+
  stat_compare_means(label = "p.format",method = "t.test",size = 4,comparisons = my_comparisons)+#ggtitle("Immune Cell Projection in Multiple Sclerosis ")+
  scale_fill_manual(values=c("#99CC66", "#FFCC66", "#E57373"))+scale_y_continuous(breaks = seq(0, 100, by = 20))+
  ylab("PD-L1 CPS")+xlab("")+
  theme( legend.position="none",legend.title = element_text( size=15, 
                                                             face="bold"),
         legend.text = element_text(size=15, 
                                    face="bold"),
         axis.text.x = element_text(size = 10, face = "bold" ),
         axis.text.y = element_text(size = 15, face = "bold" ),
         axis.title =  element_text( size = 15 , face="bold"),
         plot.title = element_text(size=20,hjust = 0.5))  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
