library(data.table)
library(mgcv)
library(mgcv.taps)
library(dplyr)
library(glue)
library(mgcViz)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
options(bitmapType = "cairo")
UKBBWithdraw <- fread("../UKBBWithdraw.csv", header = F)
UKBBKinship <- readRDS("../kinshipexclude.rds")
w1 <- readRDS("../inverse_rank_w.rds")
AgeData=w1[which(w1$Age>=60&w1$Age<=70),]
NAM_Phenotype=c("ALB","ALP",'ALT',"APOA","APOB","AST","Baso","TBL","BMI","BUN","CA","CRP","CysC","DBP",
              "Eosino","FPG","GGT","HB","HT","HBA1C","HDL","LDL","LYM","Mono","Neutro","IGF1",
              "PLA","RBC","SBP","SHBG","RET","TCh","TG","TT","TP","SCR","UA","VTD","WBC","Beer","Computer","CookedVeg","Driving","DryFruit","FEV1","FVC","PEF","Fortified","TV","ModActivity","RawVeg","FreshFruit","RedWine","WhiteWine","Bread","Tee","SleepTime","Walk","Spirits","VigActivity","Water")%>%sort(.)
NAM_No_PRS=c("Fortified","Walk","ModActivity","VigActivity","Spirits","Water")%>%sort(.)
NAM_PRS=setdiff(NAM_Phenotype,NAM_No_PRS)%>%sort(.)

phenotype_fullnames_no_prs <- c(
"Fortified"  = "Fortified Food Intake",
"ModActivity"= "Moderate Physical Activity",
"Spirits"    = "Spirits Intake",
"VigActivity"= "Vigorous Physical Activity",
"Walk"       = "Walking Time",
"Water"     = "Water Intake"
)

phenotype_fullnames_prs <- c(
"ALB"       = "Albumin",
"ALP"       = "Alkaline Phosphatase",
"ALT"       = "Alanine Aminotransferase",
"APOA"      = "Apolipoprotein A",
"APOB"      = "Apolipoprotein B",
"AST"       = "Aspartate Aminotransferase",
"Baso"      = "Basophil Count",
"Beer"      = "Beer Intake",
"BMI"       = "Body Mass Index",
"Bread"     = "Bread Intake",
"BUN"       = "Blood Urea Nitrogen",
"CA"        = "Calcium",
"Computer"  = "Computer Time",
"CookedVeg" = "Cooked Vegetable Intake",
"CRP"       = "C-reactive Protein",
"CysC"      = "Cystatin C",
"DBP"       = "Diastolic Blood Pressure",
"Driving"   = "Driving Time",
"DryFruit"  = "Dried Fruit Intake",
"Eosino"    = "Eosinophil Count",
"FEV1"      = "Forced Expiratory Volume in 1s",
"FPG"       = "Fasting Plasma Glucose",
"FreshFruit"= "Fresh Fruit Intake",
"FVC"       = "Forced Vital Capacity",
"GGT"       = "Gamma-Glutamyl Transferase",
"HB"        = "Hemoglobin",
"HBA1C"     = "Hemoglobin A1c",
"HDL"       = "HDL Cholesterol",
"HT"        = "Hematocrit",
"IGF1"      = "Insulin-like Growth Factor 1",
"LDL"       = "LDL Cholesterol",
"LYM"       = "Lymphocyte Count",
"Mono"      = "Monocyte Count",
"Neutro"    = "Neutrophil Count",
"PEF"       = "Peak Expiratory Flow",
"PLA"       = "Platelet Count",
"RawVeg"    = "Raw Vegetable Intake",
"RBC"       = "Erythrocyte Count",
"RedWine"   = "Red Wine Intake",
"RET"   = "Reticulocyte Count",
"SBP"       = "Systolic Blood Pressure",
"SCR"       = "Serum Creatinine",
"SHBG"      = "Sex Hormone-Binding Globulin",
"SleepTime" = "Sleep Duration",
"TBL"       = "Total Bilirubin",
"TCh"       = "Total Cholesterol",
"Tee" ="Tee Intake",
"TG"        = "Triglycerides",
"TP"        = "Total Protein",
"TT"        = "Total Testosterone",
"TV"        = "Television Time",
"UA"        = "Uric Acid",
"VTD"       = "Vitamin D",
"WBC"       = "Leukocytes Count",
"WhiteWine" = "White Wine Intake"
)
par(mfrow=c(1,2))
G1=G2=list()
iii=1
for(i in 1:55){
outcome=NAM_PRS[i]
form1 <- as.formula(paste0(outcome, glue("~s(Age,bs='AMatern',xt=list(getA=linearity_discontinuity,para=65))+s({outcome}_PRS,bs='cr',k=10)+Sex+Income+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
fit1=bam(form1,data=AgeData,method="REML",discrete=F)
t1=taps_wald_test(fit1)
t2=taps_score_test(fit1)
form11 <- as.formula(paste0(outcome, glue("~linearity_discontinuity(Age,65)+s({outcome}_PRS,bs='cr',k=10)+Sex+Income+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10-1")))
fit11=bam(form11,data=AgeData,method="REML",discrete=F)
g1=data.table(outcome=phenotype_fullnames_prs[i],Wald.df=t1$smooth.df,Wald.stat=t1$smooth.chisq,Wald.p=t1$smooth.pvalue,Score.df=t2$smooth.df,Score.stat=t2$smooth.stat,Score.p=t2$smooth.pvalue,jump.est=fit11$coefficients[3],jump.se=summary(fit11)$se[3],jump.p=summary(fit11)$p.pv[3])
b <- getViz(fit1)
plot_1 <- plot(sm(b, 1)) +
  l_ciPoly(mul = 5, fill = "#60c5ba", alpha = 0.5) +
  l_rug(mapping = aes(x = x), alpha = 0.5, color = "#60c5ba") +
  l_fitLine(colour = "black", size = 1) +
  theme_get() +
  xlab("Age") + ylab("f(Age)") + ggtitle(paste0("RDD: ",phenotype_fullnames_prs[i]))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "#f4f4f4")) +
  annotation_custom(
    grob = grobTree(
      textGrob(
        label = glue("Effect.df = {round(summary(fit1)$edf[1]-3,3)}\nWald.p = {format(t1$smooth.pvalue, scientific = TRUE, digits = 2)}\nScore.p = {format(t2$smooth.pvalue, scientific = TRUE, digits = 2)}"),
        x = 0.6, y = 0.96,
        hjust = 0, vjust = 1,
        gp = gpar(fontsize = 9)
      )
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )

###############################################################################################################################################################################################################
form1 <- as.formula(paste0(outcome, glue("~s(Age,bs='AMatern',xt=list(getA=piecewise_linearity,para=65))+s({outcome}_PRS,bs='cr',k=10)+Sex+Income+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
fit1=bam(form1,data=AgeData,method="REML",discrete=F)
t1=taps_wald_test(fit1)
t2=taps_score_test(fit1)
form11 <- as.formula(paste0(outcome, glue("~piecewise_linearity(Age,65)+s({outcome}_PRS,bs='cr',k=10)+Sex+Income+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10-1")))
fit11=bam(form11,data=AgeData,method="REML",discrete=F)
g2=data.table(outcome=phenotype_fullnames_prs[i],Wald.df=t1$smooth.df,Wald.stat=t1$smooth.chisq,Wald.p=t1$smooth.pvalue,Score.df=t2$smooth.df,Score.stat=t2$smooth.stat,Score.p=t2$smooth.pvalue,slopechange.est=fit11$coefficients[3],slopechange.se=summary(fit11)$se[3],slopechange.p=summary(fit11)$p.pv[3])

b <- getViz(fit1)
plot_2 <- plot(sm(b, 1)) +
  l_ciPoly(mul = 5, fill = "#f8ca00", alpha = 0.4) +
  l_rug(mapping = aes(x = x), alpha = 0.4, color = "#f8ca00") +
  l_fitLine(colour = "black", size = 1) +
  theme_get() +
  xlab("Age") + ylab("f(Age)") + ggtitle(paste0("RKD: ",phenotype_fullnames_prs[i]))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "#f4f4f4")) +
  annotation_custom(
    grob = grobTree(
      textGrob(
        label = glue("Effect.df = {round(summary(fit1)$edf[1]-2,3)}\nWald.p = {format(t1$smooth.pvalue, scientific = TRUE, digits = 2)}\nScore.p = {format(t2$smooth.pvalue, scientific = TRUE, digits = 2)}"),
        x = 0.6, y = 0.96,
        hjust = 0, vjust = 1,
        gp = gpar(fontsize = 9)
      )
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )

###############################################################################################################################################################################################################
form1 <- as.formula(paste0(outcome, glue("~s(Age,bs='cr')+s({outcome}_PRS,bs='cr',k=10)+Sex+Income+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
fit1=bam(form1,data=w1,method="REML",discrete=F)
b <- getViz(fit1)
plot_3 <- plot(sm(b, 1)) +
  l_ciPoly(mul = 5, fill = "#e85a71", alpha = 0.4) +
  l_rug(mapping = aes(x = x), alpha = 0.4, color = "#e85a71") +
  l_fitLine(colour = "black", size = 1) +
  theme_get() +
  xlab("Age") + ylab("f(Age)") + ggtitle(paste0("All Ages: ",phenotype_fullnames_prs[i]))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "#f4f4f4"))
# Save plot
png(glue("../figures/{outcome}_plot.png"), width = 12, height = 4.5, res = 600, unit = "in")
gridPrint(plot_1,plot_2,plot_3,nrow=1)
dev.off()

print(g1)
print(g2)
G1[[iii]]=g1
G2[[iii]]=g2
print(outcome)
iii=iii+1
}

for(i in 1:6){
outcome=NAM_No_PRS[i]
form1 <- as.formula(paste0(outcome, glue("~s(Age,bs='AMatern',xt=list(getA=linearity_discontinuity,para=65))+Sex+Income+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
fit1=bam(form1,data=AgeData,method="REML")
t1=taps_wald_test(fit1)
t2=taps_score_test(fit1)
form11 <- as.formula(paste0(outcome, glue("~linearity_discontinuity(Age,65)+Sex+Income+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10-1")))
fit11=bam(form11,data=AgeData,method="REML")
g1=data.table(outcome=phenotype_fullnames_no_prs[i],Wald.df=t1$smooth.df,Wald.stat=t1$smooth.chisq,Wald.p=t1$smooth.pvalue,Score.df=t2$smooth.df,Score.stat=t2$smooth.stat,Score.p=t2$smooth.pvalue,jump.est=fit11$coefficients[3],jump.se=summary(fit11)$se[3],jump.p=summary(fit11)$p.pv[3])
b <- getViz(fit1)
plot_1 <- plot(sm(b, 1)) +
  l_ciPoly(mul = 5, fill = "#60c5ba", alpha = 0.5) +
  l_rug(mapping = aes(x = x), alpha = 0.5, color = "#60c5ba") +
  l_fitLine(colour = "black", size = 1) +
  theme_get() +
  xlab("Age") + ylab("f(Age)") + ggtitle(paste0("RDD: ",phenotype_fullnames_no_prs[i]))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "#f4f4f4")) +
  annotation_custom(
    grob = grobTree(
      textGrob(
        label = glue("Effect.df = {round(summary(fit1)$edf[1]-3,3)}\nWald.p = {format(t1$smooth.pvalue, scientific = TRUE, digits = 2)}\nScore.p = {format(t2$smooth.pvalue, scientific = TRUE, digits = 2)}"),
        x = 0.6, y = 0.96,
        hjust = 0, vjust = 1,
        gp = gpar(fontsize = 9)
      )
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )
###############################################################################################################################################################################################################
form1 <- as.formula(paste0(outcome, glue("~s(Age,bs='AMatern',xt=list(getA=piecewise_linearity,para=65))+Sex+Income+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
fit1=bam(form1,data=AgeData,method="REML")
t1=taps_wald_test(fit1)
t2=taps_score_test(fit1)
form11 <- as.formula(paste0(outcome, glue("~piecewise_linearity(Age,65)+Sex+Income+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10-1")))
fit11=bam(form11,data=AgeData,method="REML")
g2=data.table(outcome=phenotype_fullnames_no_prs[i],Wald.df=t1$smooth.df,Wald.stat=t1$smooth.chisq,Wald.p=t1$smooth.pvalue,Score.df=t2$smooth.df,Score.stat=t2$smooth.stat,Score.p=t2$smooth.pvalue,slopechange.est=fit11$coefficients[3],slopechange.se=summary(fit11)$se[3],slopechange.p=summary(fit11)$p.pv[3])
b <- getViz(fit1)
plot_2 <- plot(sm(b, 1)) +
  l_ciPoly(mul = 5, fill = "#f8ca00", alpha = 0.4) +
  l_rug(mapping = aes(x = x), alpha = 0.4, color = "#f8ca00") +
  l_fitLine(colour = "black", size = 1) +
  theme_get() +
  xlab("Age") + ylab("f(Age)") + ggtitle(paste0("RKD: ",phenotype_fullnames_no_prs[i]))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "#f4f4f4")) +
  annotation_custom(
    grob = grobTree(
      textGrob(
        label = glue("Effect.df = {round(summary(fit1)$edf[1]-2,3)}\nWald.p = {format(t1$smooth.pvalue, scientific = TRUE, digits = 2)}\nScore.p = {format(t2$smooth.pvalue, scientific = TRUE, digits = 2)}"),
        x = 0.6, y = 0.96,
        hjust = 0, vjust = 1,
        gp = gpar(fontsize = 9)
      )
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )

###############################################################################################################################################################################################################
form1 <- as.formula(paste0(outcome, glue("~s(Age,bs='cr')+Sex+Income+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
fit1=bam(form1,data=w1,method="REML",discrete=F)
b <- getViz(fit1)
plot_3 <- plot(sm(b, 1)) +
  l_ciPoly(mul = 5, fill = "#e85a71", alpha = 0.4) +
  l_rug(mapping = aes(x = x), alpha = 0.4, color = "#e85a71") +
  l_fitLine(colour = "black", size = 1) +
  theme_get() +
  xlab("Age") + ylab("f(Age)") + ggtitle(paste0("All Ages: ",phenotype_fullnames_prs[i]))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "#f4f4f4"))

png(glue("../figures/{outcome}_cbind_plot.png"), width = 12, height = 4.5, res = 600, unit = "in")
gridPrint(plot_1,plot_2,plot_3,nrow=1)
dev.off()
print(g1)
print(g2)
G1[[iii]]=g1
G2[[iii]]=g2
print(outcome)
iii=iii+1
}

G1=do.call(rbind,G1)
G2=do.call(rbind,G2)


