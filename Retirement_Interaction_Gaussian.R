library(data.table)
library(dplyr)
library(glue)
library(mgcv)
library(mgcv.taps)
library(ggplot2)
library(glmnet)
library(rcompanion)
library(mgcViz)
library(grid)
library(scales)
library(gratia)
options(bitmapType="cairo")
inverse_rank_normalize <- function(x, ties.method = "average") {
not_na <- !is.na(x)
n <- sum(not_na)
ranks <- rep(NA_real_, length(x))
ranks[not_na] <- rank(x[not_na], ties.method = ties.method)
quantiles <- rep(NA_real_, length(x))
quantiles[not_na] <- (ranks[not_na] - 0.5) / n
qnorm(quantiles)
}

cutoff_quantile=function(w,d=0.01){
ind1=which(w<quantile(w,d))
ind2=which(w>quantile(w,1-d))
ind=union(ind1,ind2)
w[ind]=NA
return(w)
}


UKBBWithdraw <- fread("~/1000G/UKBBWithdraw.csv", header = F)
UKBBKinship <- readRDS("~/MRJones_CRE/kinshipexclude.rds")
w1 <- readRDS("~/mgcv_maps_simulation/PRS_UKBB.rds") %>% as.data.frame(.)
colnames(w1)[-1]=paste0(colnames(w1)[-1],"_PRS")
w2 <- readRDS("~/mgcv_maps_simulation/UKBB_Phenotype_Data.rds")
NAM_Phenotype=colnames((w2))
w1 <- merge(w2, w1, by = "ID")
w2 <- readRDS("~/mgcv_maps_simulation/UKBB_Disease.rds")
NAM_Disease=colnames((w2))
w1 <- merge(w2, w1, by = "ID")
w1$SBP=ifelse(w1$BPDrug==1,w1$SBP+15,w1$SBP)
w1$DBP=ifelse(w1$BPDrug==1,w1$DBP+10,w1$DBP)
w1$PP=w1$SBP-w1$DBP
w2=readRDS("~/MRJones_CRE/PRS/NewBP_PRS_50000.rds")%>%dplyr::select(ID,PP_PRS=PP_All)
w1=merge(w1,w2,by="ID")
ind <- which(w1$ID %in% union(UKBBWithdraw$V1, UKBBKinship))
w1 <- w1[-ind, ]
eurind <- readRDS("~/Mr.Jones/indeur.rds")
w1 <- w1[which(w1$ID %in% eurind), ]

w1=w1%>%dplyr::filter(Age>63&Age<67)

NAM_Disease=NAM_Disease[c(2,3,6,8,9,10,11,12)]
NAM_Phenotype=c("ALB","ALP",'ALT',"APOA","APOB","AST","TBL","BMI","BUN","CA","CRP","CysC","DBP","FPG","GGT","HB","HT","HBA1C","HDL","HEG","LDL","LYM","Mono","Neutro","IGF1","PLA","PP","RBC","SBP","SHBG","RET","TCh","TG","TT","TP","UA","VTD","WBC")
phenotype_fullnames <- c(
ALB   = "Albumin",
ALP   = "Alkaline Phosphatase",
ALT   = "Alanine Aminotransferase",
APOA  = "Apolipoprotein A",
APOB  = "Apolipoprotein B",
AST   = "Aspartate Aminotransferase",
TBL   = "Total Bilirubin",
BMI   = "Body Mass Index",
BUN   = "Blood Urea Nitrogen",
CA    = "Calcium",
CRP   = "C-reactive Protein",
CysC  = "Cystatin C",
DBP   = "Diastolic Blood Pressure",
FPG   = "Fasting Plasma Glucose",
GGT   = "Gamma-Glutamyl Transferase",
HB    = "Hemoglobin",
HT    = "Hematocrit",
HBA1C = "Hemoglobin A1c",
HDL = "HDL Cholesterol",
HEG   = "Height",
LDL = "LDL Cholesterol",
LYM   = "Lymphocyte Count",
Mono  = "Monocyte Count",
Neutro= "Neutrophil Count",
IGF1  = "Insulin-like Growth Factor 1",
PLA   = "Platelet Count",
PP   = "Pulse Pressure",
RBC   = "Erythrocyte Count",
SBP   = "Systolic Blood Pressure",
SHBG  = "Sex Hormone-Binding Globulin",
RET   = "Reticulocyte Count",
TCh = "Total Cholesterol",
TG = "Triglyceride",
TT    = "Total Testosterone",
TP    = "Total Protein",
UA    = "Uric Acid",
VTD   = "Vitamin D",
WBC   = "Leukocytes Count"
)

G1=G2=list()
for (i in 1:38) {
trait <- NAM_Phenotype[i]
if(trait!="LDL"|trait!="TCh"){
wi=w1%>%dplyr::select(trait,paste0(trait,"_PRS"),Age,Sex,paste0("PC",1:10))%>%as.data.frame(.)%>%na.omit(.)
names(wi)[1]="Outcome"
wi$Outcome=inverse_rank_normalize(wi$Outcome);
formula_Outcome <- as.formula(glue("Outcome ~ s(Age,by={paste0(trait,'_PRS')},bs='AMatern',k=8,xt=list(getA=linearity_discontinuity,para=65))+s(Age,bs='cr')+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"))
}else{
wi=w1%>%dplyr::select(trait,paste0(trait,"_PRS"),Age,Sex,CholDrug,paste0("PC",1:10))%>%as.data.frame(.)%>%na.omit(.)
names(wi)[1]="Outcome"
wi$Outcome=inverse_rank_normalize(wi$Outcome);
formula_Outcome <- as.formula(glue("Outcome ~ s(Age,by={paste0(trait,'_PRS')},bs='AMatern',k=8,xt=list(getA=linearity_discontinuity,para=65))+s(Age,bs='cr')+Sex+CholDrug+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"))
}
fit_Outcome=bam(formula_Outcome,method="REML",discrete=F,data=wi)
wald_Outcome=taps_wald_test(fit_Outcome)
score_Outcome=taps_score_test(fit_Outcome)
g_Outcome=data.frame(outcome=trait,edf=wald_Outcome$smooth.df,int.est=coef(fit_Outcome)[glue("s(Age):{paste0(trait,'_PRS')}.3")],int.se=summary(fit_Outcome)$se[glue("s(Age):{paste0(trait,'_PRS')}.3")],wald.p=wald_Outcome$smooth.pvalue,score.p=score_Outcome$smooth.pvalue)

b <- getViz(fit_Outcome)
plot_1 <- plot(sm(b, 1)) +
l_ciPoly(mul = 5, fill = "#60c5ba", alpha = 0.35) +
l_rug(mapping = aes(x = x), alpha = 0.35, color = "#60c5ba") +
l_fitLine(colour = "black", size = 1) +
geom_hline(yintercept=0,linetype="dashed")+
theme_get() + ggtitle(glue("Outcome: {trait}\n Jump at Age = 65"))+
xlab("Age") + ylab((glue('{trait} PRS x s(Age)'))) +
theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
panel.background = element_rect(fill = "#f4f4f4")) +
annotation_custom(
grob = grobTree(
textGrob(
label = glue("Effect.df = {round(summary(fit_Outcome)$edf[1],3)}\nWald.p = {format(g_Outcome$wald.p, scientific = TRUE, digits = 2)}\nScore.p = {format(g_Outcome$score.p, scientific = TRUE, digits = 2)}"),
x = 0.8, y = 0.3,
hjust = 0, vjust = 1,
gp = gpar(fontsize = 11)
)
),
xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
)

G1[[i]]=g_Outcome
png(glue("~/mgcv_taps_varying/figures/{trait}_RDD_quan.png"), width = 7, height = 7, res = 600, unit = "in")
gridPrint(plot_1)
dev.off()
#############################################################################################################
if(trait!="LDL"|trait!="TCh"){
wi=w1%>%dplyr::select(trait,paste0(trait,"_PRS"),Age,Sex,paste0("PC",1:10))%>%as.data.frame(.)%>%na.omit(.)
names(wi)[1]="Outcome"
wi$Outcome=inverse_rank_normalize(wi$Outcome);
formula_Outcome <- as.formula(glue("Outcome ~ s(Age,by={paste0(trait,'_PRS')},bs='AMatern',k=8,xt=list(getA=piecewise_linearity,para=65))+s(Age,bs='cr')+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"))
}else{
wi=w1%>%dplyr::select(trait,paste0(trait,"_PRS"),Age,Sex,CholDrug,paste0("PC",1:10))%>%as.data.frame(.)%>%na.omit(.)
names(wi)[1]="Outcome"
wi$Outcome=inverse_rank_normalize(wi$Outcome);
formula_Outcome <- as.formula(glue("Outcome ~ s(Age,by={paste0(trait,'_PRS')},bs='AMatern',k=8,xt=list(getA=piecewise_linearity,para=65))+s(Age,bs='cr')+Sex+CholDrug+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"))
}
fit_Outcome=bam(formula_Outcome,method="REML",discrete=F,data=wi)
wald_Outcome=taps_wald_test(fit_Outcome)
score_Outcome=taps_score_test(fit_Outcome)
g_Outcome=data.frame(outcome=trait,edf=wald_Outcome$smooth.df,int.est=coef(fit_Outcome)[glue("s(Age):{paste0(trait,'_PRS')}.3")],int.se=summary(fit_Outcome)$se[glue("s(Age):{paste0(trait,'_PRS')}.3")],wald.p=wald_Outcome$smooth.pvalue,score.p=score_Outcome$smooth.pvalue)

b <- getViz(fit_Outcome)
plot_1 <- plot(sm(b, 1)) +
l_ciPoly(mul = 5, fill = "#fcbe32", alpha = 0.35) +
l_rug(mapping = aes(x = x), alpha = 0.35, color = "#fcbe32") +
l_fitLine(colour = "black", size = 1) +
geom_hline(yintercept=0,linetype="dashed")+
theme_get() + ggtitle(glue("Outcome: {trait}\n Slope Change at Age = 65"))+
xlab("Age") + ylab((glue('{trait} PRS x s(Age)'))) +
theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
panel.background = element_rect(fill = "#f4f4f4")) +
annotation_custom(
grob = grobTree(
textGrob(
label = glue("Effect.df = {round(summary(fit_Outcome)$edf[1],3)}\nWald.p = {format(g_Outcome$wald.p, scientific = TRUE, digits = 2)}\nScore.p = {format(g_Outcome$score.p, scientific = TRUE, digits = 2)}"),
x = 0.8, y = 0.3,
hjust = 0, vjust = 1,
gp = gpar(fontsize = 11)
)
),
xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
)

G2[[i]]=g_Outcome
png(glue("~/mgcv_taps_varying/figures/{trait}_RKD_quan.png"), width = 7, height = 7, res = 600, unit = "in")
gridPrint(plot_1)
dev.off()

}
G1=do.call(rbind,G1)
G2=do.call(rbind,G2)
write.csv(G1,"~/mgcv_taps_varying/Gaussian_Exposure_RDD.csv")
write.csv(G2,"~/mgcv_taps_varying/Gaussian_Exposure_RKD.csv")

