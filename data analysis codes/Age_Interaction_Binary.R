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
UKBBWithdraw <- fread("~/1000G/UKBBWithdraw.csv", header = F)
UKBBKinship <- readRDS("~/MRJones_CRE/kinshipexclude.rds")
w1 <- readRDS("~/mgcv_maps_simulation/PRS_UKBB.rds") %>% as.data.frame(.)
colnames(w1)[-1]=paste0(colnames(w1)[-1],"_PRS")
w2 <- readRDS("~/mgcv_maps_simulation/UKBB_Phenotype_Data.rds")
NAM_Phenotype=colnames((w2))
w1 <- merge(w2, w1, by = "ID")
w2 <- readRDS("~/mgcv_maps_simulation/UKBB_Disease.rds")%>%dplyr::select(ID,CAD,Stroke,HF,AF,PAD,T2D,SMK,paste0("PC",1:10))
NAM_Disease=colnames((w2))
w1 <- merge(w2, w1, by = "ID")
w2 = readRDS("~/SuSiE4X_Simulation/Drinking.rds")
w2$DRNK=w2$Daily+w2$Frequently
w1 <- merge(w2, w1, by = "ID")
w2=readRDS("~/mgcv_maps_simulation/SMK_DRNK.rds")%>%dplyr::select(ID,DRNK_PRS=DRNK)
w1 <- merge(w2, w1, by = "ID")

ind <- which(w1$ID %in% union(UKBBWithdraw$V1, UKBBKinship))
w1 <- w1[-ind, ]
eurind <- readRDS("~/Mr.Jones/indeur.rds")
w1 <- w1[which(w1$ID %in% eurind), ]
NAM_Disease=NAM_Disease[2:7]
NAM_Disease=c(NAM_Disease,"DRNK")
NAM_Phenotype=c("ALB","ALP",'ALT',"APOA","APOB","AST","TBL","BMI","BUN","CA","CRP","CysC","DBP","FPG","GGT","HB","HT","HBA1C","HDL","HEG","LDL","LYM","Mono","Neutro","IGF1","PLA","PP","RBC","SBP","SHBG","RET","TCh","TG","TT","TP","UA","VTD","WBC")
disease_fullnames <- c(
CAD          = "Coronary Artery Disease",
Stroke       = "Ischemic Stroke",
HF           = "Heart Failure",
AF           = "Atrial Fibrillation",
PAD          = "Peripheral Artery Disease",
T2D          = "Type 2 Diabetes",
SMK          = "Smoking Initiation",
DRNK         = "Frequently Drink"
)
###############################################################################
test_results=G<- list()

###############################################################################

for (i in 10:10) {
trait <- NAM_Disease[i]

wi=w1%>%dplyr::select(trait,paste0(trait,"_PRS"),Age,Sex,paste0("PC",1:10))%>%as.data.frame(.)
formula <- as.formula(glue("{trait} ~ s(Age,by={paste0(trait,'_PRS')},bs='AMatern',k=10) + s(Age,bs='cr') + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
fit <- bam(formula, data = w1, method = "REML",discrete = F, family = "binomial")

# Store Wald test result with trait name
t1 <- taps_wald_test(fit)[,6:8] %>% dplyr::select(df=smooth.df,stat=smooth.chisq,p=smooth.pvalue)
t2 <- taps_score_test(fit)[,2:4] %>% dplyr::select(df=smooth.df,stat=smooth.stat,p=smooth.pvalue)
t1$Test="Wald Test"
t2$Test="Score Test"
t=rbind(t1,t2)
t$Trait=disease_fullnames[i]
test_results[[i]]=t

formula1 <- as.formula(glue("{trait} ~ {paste0(trait,'_PRS')} + s(Age,bs='cr') + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
fit1 <- bam(formula1, data = w1, method = "REML",discrete = F, family = "binomial")

test=anova(fit1,fit)
g_Outcome=data.frame(outcome=trait,edf=t1$df,int.est=coef(fit)[glue("s(Age):{paste0(trait,'_PRS')}.2")],int.se=summary(fit)$se[glue("s(Age):{paste0(trait,'_PRS')}.2")],wald.p=t1$p,score.p=t2$p)
G[[i]]=g_Outcome

b <- getViz(fit)
plot_obj <- plot(sm(b, 1)) +
l_ciPoly(mul = 5, fill = "#60c5ba", alpha = 0.25) +
l_rug(mapping = aes(x = x), alpha = 0.25, color = "#60c5ba") +
l_fitLine(colour = "black", size = 1) +
theme_get() +
theme_get() + ggtitle(glue("Outcome: {trait}\nVarying Coefficient: Age"))+
xlab("Age") + ylab((glue('{trait} PRS x s(Age)'))) +
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
png(glue("~/mgcv_taps_varying/figures/{trait}_plot.png"), width = 7, height = 7, res = 600, unit = "in")
print(plot_obj)
dev.off()

test_Outcome=anova(fit1,fit)
g_Outcome=data.frame(Trait=trait,v1=summary(fit1)$r.sq,v2=summary(fit)$r.sq,Fstat=test_Outcome$`Pr(>Chi)`[2],pv=test_Outcome$`Pr(>Chi)`[2])
g_Outcome$v1 <- percent(g_Outcome$v1, accuracy = 0.001)
g_Outcome$v2 <- percent(g_Outcome$v2, accuracy = 0.001)
test_results[[i]]=g_Outcome

print(trait)
}

G=do.call(rbind,G)
write.csv(G,"~/mgcv_taps_varying/Binary_Exposure.csv")
test_results=do.call(rbind,test_results)
saveRDS(test_results, "~/mgcv_taps_varying/Binary_Test_Results.rds")
