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

get_quantile_indices <- function(x, n = 1000) {
x <- as.numeric(x)
probs <- seq(0, 1, length.out = n)
q_vals <- quantile(x, probs = probs, names = FALSE)
idx <- sapply(q_vals, function(q) {
which.min(abs(x - q))
})
return(unique(idx))
}

# Data preparation (unchanged)
UKBBWithdraw <- fread("../UKBBWithdraw.csv", header = F)
UKBBKinship <- readRDS("../kinshipexclude.rds")
w1 <- readRDS("../PRS_UKBB.rds") %>% as.data.frame(.)
colnames(w1)[-1]=paste0(colnames(w1)[-1],"_PRS")
w2 <- readRDS("../UKBB_Phenotype_Data.rds")
NAM_Phenotype=colnames((w2))
w1 <- merge(w2, w1, by = "ID")
w2 <- readRDS("../UKBB_Disease.rds")
NAM_Disease=colnames((w2))
w1 <- merge(w2, w1, by = "ID")
ind <- which(w1$ID %in% union(UKBBWithdraw$V1, UKBBKinship))
w1 <- w1[-ind, ]
eurind <- readRDS("../indeur.rds")
w1 <- w1[which(w1$ID %in% eurind), ]
names(w1)[61]="BUN"
NAM_Disease=NAM_Disease[c(2,3,6,8,9,10,11,12)]
NAM_Phenotype=c("ALB","ALP",'ALT',"APOA","APOB","AST","TBL","BMI","BUN","CA","CRP","CysC","DBP",
             "FPG","GGT","HB","HT","HBA1C","HDL","HEG","LDL","LYM","Mono","Neutro","IGF1",
             "PLA","RBC","SBP","SHBG","RET","TCh","TG","TT","TP","UA","VTD","WBC","SCR","Eosino","Baso")
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
HDL   = "HDL Cholesterol",
HEG   = "Height",
LDL   = "LDL Cholesterol",
LYM   = "Lymphocyte Count",
Mono  = "Monocyte Count",
Neutro= "Neutrophil Count",
IGF1  = "Insulin-like Growth Factor 1",
PLA   = "Platelet Count",
RBC   = "Erythrocyte Count",
SBP   = "Systolic Blood Pressure",
SHBG  = "Sex Hormone-Binding Globulin",
RET   = "Reticulocyte Count",
TCh   = "Total Cholesterol",
TG    = "Triglycerides",
TT    = "Total Testosterone",
TP    = "Total Protein",
UA    = "Uric Acid",
VTD   = "Vitamin D",
WBC   = "Leukocytes Count",
SCR   = "Serum Creatinine",
Eosino= "Eosinophil Count",
Baso  = "Basophil Count"
)

test_results1=test_results2= test_results3= test_results4<- list()
G1=G2=G3=G4=list()
###############################################################################

for (i in 1:37) {
trait <- NAM_Phenotype[i]
wi <- w1 %>%
  dplyr::select(
    trait,
    paste0(trait, "_PRS"),
    paste0(trait, "_Raw_PRS"),
    Age,
    Sex,
    paste0("PC", 1:10)
  ) %>%
  as.data.frame()%>%na.omit(.)
wi$quan=inverse_rank_normalize(wi[,1])
wi[,1]=cutoff_quantile(wi[,1])
formula <- as.formula(glue("{trait} ~ s({paste0(trait,'_Raw_PRS')},bs='AMatern',k=10) + s(Age,bs='cr') + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
fit <- bam(formula, data = wi, method = "REML",discrete = F)

# Store Wald test result with trait name
t1 <- taps_wald_test(fit)[,6:8] %>% dplyr::select(df=smooth.df,stat=smooth.chisq,p=smooth.pvalue)
t2 <- taps_score_test(fit)[,2:4] %>% dplyr::select(df=smooth.df,stat=smooth.stat,p=smooth.pvalue)
t1$Test="Wald Test"
t2$Test="Score Test"
t=rbind(t1,t2)
t$Trait=phenotype_fullnames[i]
test_results1[[i]]=t

formula1 <- as.formula(glue("{trait} ~ {paste0(trait,'_Raw_PRS')} + s(Age,bs='cr') + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
fit1 <- bam(formula1, data = wi, method = "REML",discrete = F)
test=anova(fit1,fit)
g=data.frame(trait=phenotype_fullnames[i],v1=summary(fit1)$r.sq,v2=summary(fit)$r.sq,Fstat=test$F[2],pv=test$`Pr(>F)`[2])
g$v1 <- percent(g$v1, accuracy = 0.01)
g$v2 <- percent(g$v2, accuracy = 0.01)
G1[[i]]=g
# Generate plot stfit# Generate plot styled like plot2
b <- getViz(fit)
plot_1 <- plot(sm(b, 1)) +
l_ciPoly(mul = 5, fill = "#60c5ba", alpha = 0.25) +
l_rug(mapping = aes(x = x), alpha = 0.25, color = "#60c5ba") +
l_fitLine(colour = "black", size = 1) +
theme_get() + ggtitle("Raw Outcome")+
xlab(paste0(phenotype_fullnames[i]," PRS",glue(" ({trait} PRS)"))) + ylab(glue('f({trait} PRS)')) +
theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "#f4f4f4")) +
annotation_custom(
grob = grobTree(
  textGrob(
    label = glue("Effect.df = {round(summary(fit)$edf[1]-1,3)}\nWald.p = {format(t1$p, scientific = TRUE, digits = 2)}\nScore.p = {format(t2$p, scientific = TRUE, digits = 2)}"),
    x = 0.03, y = 0.965,
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 11)
  )
),
xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
)

############################################################################################################
formula <- as.formula(glue("quan ~ s({paste0(trait,'_PRS')},bs='AMatern',k=10) + s(Age,bs='cr') + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
fit <- bam(formula, data = wi, method = "REML",discrete = F)

# Store Wald test result with trait name
t1 <- taps_wald_test(fit)[,6:8] %>% dplyr::select(df=smooth.df,stat=smooth.chisq,p=smooth.pvalue)
t2 <- taps_score_test(fit)[,2:4] %>% dplyr::select(df=smooth.df,stat=smooth.stat,p=smooth.pvalue)
t1$Test="Wald Test"
t2$Test="Score Test"
t=rbind(t1,t2)
t$Trait=phenotype_fullnames[i]
test_results2[[i]]=t

formula1 <- as.formula(glue("quan ~ {paste0(trait,'_PRS')} + s(Age,bs='cr') + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
fit1 <- bam(formula1, data = wi, method = "REML",discrete = F)
test=anova(fit1,fit)
g=data.frame(trait=phenotype_fullnames[i],v1=summary(fit1)$r.sq,v2=summary(fit)$r.sq,Fstat=test$F[2],pv=test$`Pr(>F)`[2])
g$v1 <- percent(g$v1, accuracy = 0.01)
g$v2 <- percent(g$v2, accuracy = 0.01)
G2[[i]]=g
# Generate plot stfit# Generate plot styled like plot2
b <- getViz(fit)
plot_2 <- plot(sm(b, 1)) +
  l_ciPoly(mul = 5, fill = "#f8ca00", alpha = 0.25) +
  l_rug(mapping = aes(x = x), alpha = 0.25, color = "#f8ca00") +
  l_fitLine(colour = "black", size = 1) +
  theme_get() + ggtitle("Inverse-Rank Normalized Outcome")+
  xlab(paste0(phenotype_fullnames[i]," PRS",glue(" ({trait} PRS)"))) + ylab(glue('f({trait} PRS)')) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "#f4f4f4")) +
  annotation_custom(
    grob = grobTree(
      textGrob(
        label = glue("Effect.df = {round(summary(fit)$edf[1]-1,3)}\nWald.p = {format(t1$p, scientific = TRUE, digits = 2)}\nScore.p = {format(t2$p, scientific = TRUE, digits = 2)}"),
        x = 0.03, y = 0.965,
        hjust = 0, vjust = 1,
        gp = gpar(fontsize = 11)
      )
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )

# Save plot
png(glue("../figures/{trait}_plot.png"), width = 8, height = 4, res = 600, unit = "in")
gridPrint(plot_1,plot_2,nrow=1)
dev.off()

png(glue("../figures/{trait}_Raw_plot.png"), width = 4, height = 4, res = 600, unit = "in")
plot_1
dev.off()
print(trait)
}

G1=do.call(rbind,G1)
G2=do.call(rbind,G2)
test_results1=do.call(rbind,test_results1)
test_results2=do.call(rbind,test_results2)

