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
options(bitmapType="cairo")
inverse_rank_normalize <- function(x, ties.method = "average") {
ranks <- rank(x, ties.method = ties.method)
quantiles <- (ranks - 0.5) / length(x)
qnorm(quantiles)
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
NAM_Disease=NAM_Disease[c(2,3,7,8,9,10,11,12)]
NAM_Phenotype=c("ALB","ALP",'ALT',"APOA","APOB","AST","Baso","TBL","BMI","BUN","CA","CRP","CysC","DBP",
              "Eosino","FPG","GGT","HB","HT","HBA1C","HDL","HEG","LDL","LYM","Mono","Neutro","IGF1",
              "PLA","RBC","SBP","SHBG","RET","TCh","TG","TT","TP","SCR","UA","VTD","WBC")
disease_fullnames <- c(
Hypertension = "Hypertension",
CAD          = "Coronary Artery Disease",
Stroke       = "Ischemic Stroke",
HF           = "Heart Failure",
AF           = "Atrial Fibrillation",
PAD          = "Peripheral Artery Disease",
T2D          = "Type 2 Diabetes",
SMK          = "Smoking Initiation"
)
#i=which(names(w1)=="Stroke_PRS")
#names(w1)[i]="CVD_PRS"
test_results <- list()
G=list()
###############################################################################

for (i in 1:8) {
trait <- NAM_Disease[i]
if(i==1){
wi=w1%>%dplyr::select(trait,SBP_PRS,DBP_PRS,Age,Sex,paste0("PC",1:10))%>%as.data.frame(.)
}else{
wi=w1%>%dplyr::select(trait,paste0(trait,"_PRS"),Age,Sex,paste0("PC",1:10))%>%as.data.frame(.)
}

if(i!=1){
formula <- as.formula(glue("{trait} ~ s({paste0(trait,'_PRS')},bs='AMatern',k=10) + s(Age,bs='cr') + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
fit <- bam(formula, data = w1, method = "REML",discrete = F, family = "binomial")
}else{
formula <- as.formula(glue("{trait} ~ s(SBP_PRS,DBP_PRS,bs='A2Matern',k=10) + s(Age,bs='cr') + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
fit <- bam(formula, data = w1, method = "REML",discrete = F, family = "binomial")
}
# Store Wald test result with trait name
t1 <- taps_wald_test(fit)[,6:8] %>% dplyr::select(df=smooth.df,stat=smooth.chisq,p=smooth.pvalue)
t2 <- taps_score_test(fit)[,2:4] %>% dplyr::select(df=smooth.df,stat=smooth.stat,p=smooth.pvalue)
t1$Test="Wald Test"
t2$Test="Score Test"
t=rbind(t1,t2)
t$Trait=disease_fullnames[i]
test_results[[i]]=t

if(i!=1){
formula1 <- as.formula(glue("{trait} ~ {paste0(trait,'_PRS')} + s(Age,bs='cr') + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
fit1 <- bam(formula1, data = w1, method = "REML",discrete = F, family = "binomial")
}else{
formula1 <- as.formula(glue("{trait} ~ SBP_PRS*DBP_PRS + s(Age,bs='cr') + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
fit1 <- bam(formula1, data = w1, method = "REML",discrete = F, family = "binomial")
}
test=anova(fit1,fit)
g=data.frame(trait=disease_fullnames[i],v1=summary(fit1)$r.sq,v2=summary(fit)$r.sq,Chistat=test$Deviance[2],pv=test$`Pr(>Chi)`[2])
g$v1 <- percent(g$v1, accuracy = 0.01)
g$v2 <- percent(g$v2, accuracy = 0.01)
G[[i]]=g
# Generate plot stfit# Generate plot styled like plot2
b <- getViz(fit)
if(i!=1){
plot_obj <- plot(sm(b, 1)) +
  l_ciPoly(mul = 5, fill = "#60c5ba", alpha = 0.25) +
  l_rug(mapping = aes(x = x), alpha = 0.25, color = "#60c5ba") +
  l_fitLine(colour = "black", size = 1) +
  theme_get() +
  xlab(paste0(disease_fullnames[i]," PRS",glue(" ({trait} PRS)"))) + ylab(glue('f({trait} PRS)')) +
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
png(glue("../figures/{trait}_plot.png"), width = 4, height = 4, res = 600, unit = "in")
print(plot_obj)
dev.off()
}else{
myvis.gam <- readRDS("../myvis.gam.rds")
jet.colors=colorRampPalette(c("#D7FFF1", "#60c5ba", "#004e66"))
png("../figures/HT_mesh.png",width=7,height=7,res=300,unit="in")
myvis.gam(fit,
          view = c("SBP_PRS", "DBP_PRS"),
          plot.type = "persp",
          theta = 300, phi = 30,
          xlab = "SBP PRS",  # X-axis label
          ylab = "DBP PRS",  # Y-axis label
          zlab = "f(SBP PRS, DBP PRS)",  # Z-axis (response) label
          main = "T. Hypertension",  # Title
          color = "jet")
dev.off()
}
print(trait)
}
G=do.call(rbind,G)
test_results=do.call(rbind,test_results)
