
library(data.table)
library(dplyr)
library(ggplot2)
library(readr)
library(forcats)

options(bitmapType = "cairo") 


Gaussian_Test_Results <- readRDS("~/mgcv_maps_simulation/Gaussian_Test_Results.rds")%>%mutate(Type="Gaussian")
G1=Gaussian_Test_Results[which(Gaussian_Test_Results$Test=="Wald Test"),]
G2=Gaussian_Test_Results[which(Gaussian_Test_Results$Test=="Score Test"),]
G1=G1%>%dplyr::select(Trait,Type,Edf=df,`Wald Test`=p)
G1$`Score Test`=G2$p
Binary_Test_Results <- readRDS("~/mgcv_maps_simulation/Binary_Test_Results.rds")%>%mutate(Type="Binary")
B1=Binary_Test_Results[which(Binary_Test_Results$Test=="Wald Test"),]
B2=Binary_Test_Results[which(Binary_Test_Results$Test=="Score Test"),]
B1=B1%>%dplyr::select(Trait,Type,Edf=df,`Wald Test`=p)
B1$`Score Test`=B2$p

Test_Results <- do.call(rbind,list(G1,B1))%>%arrange(Type,Trait)
Test_Results = rbind(Binary_Test_Results,Gaussian_Test_Results)

trait_category <- read_csv("trait_category_mapping.csv")  
trait_category$Trait[3]="Ischemic Stroke"

score_df <- Test_Results %>%
  filter(Test == "Score Test") %>%
  mutate(
    p = ifelse(p == 0, 1e-322, p),
    log10P = -log10(p)
  ) %>%
  left_join(trait_category, by = "Trait")

category_order <- c(
  "Cardiovascular Diseases", "Diabetes Related", "Basic Traits",
  "Lipid","Blood Pressure", "Liver Function", "Renal Function",
  "Blood Cell", "Hormone Related", "Other"
)
score_df$Category <- factor(score_df$Category, levels = category_order)

score_df <- score_df %>%
  arrange(Category, desc(log10P)) %>%
  mutate(Trait = factor(Trait, levels = unique(Trait)))

brbg_colors <- c(
  "#003c30", "#543005",
  "#01665e", "#8c510a",
  "#35978f", "#bf812d",
  "#80cdc1", "#dfc27d",
  "#c7eae5", "#f6e8c3"
)

svg("~/mgcv_maps_simulation/figures/manhatten.svg",width=13,height=5.5)
ggplot(score_df, aes(x = Trait, y = log10P, color = Category)) +
  geom_segment(aes(xend = Trait, y = 0, yend = log10P), linewidth = 1.2) +
  geom_segment(aes(xend = Trait, y = 0, yend = log10P), linewidth = 0.2, color = "black") +
  geom_point(size = 2.8) +
  geom_point(size = 0.8, color = "black") +
  scale_color_manual(values = brbg_colors)+
  labs(
    x = "Trait",
    y = expression(-log[10](italic(P)))
  ) +
  geom_text(
    data = score_df[score_df$log10P > -log10(0.05/50), ],
    aes(x = Trait, y = log10P + 10, label = "*"),
    color = "black",
    size = 5,
    vjust = 0.5
  )+
  guides(color = guide_legend(title = NULL)) + 
  guides(color = guide_legend(title = NULL, nrow = 1)) +  
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),               
    legend.key.size = unit(0.5, "lines"),               
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)
  )
dev.off()

Test_Results <- readRDS("~/mgcv_maps_simulation/Gaussian_Test_Results_quantile.rds")%>%mutate(Type="Gaussian")
G1=Test_Results[which(Test_Results$Test=="Wald Test"),]
G2=Test_Results[which(Test_Results$Test=="Score Test"),]
G1=G1%>%dplyr::select(Trait,Type,Edf=df,`Wald Test`=p)
G1$`Score Test`=G2$p
G1=arrange(G1,Type,Trait)
score_df <- Test_Results %>%
  filter(Test == "Score Test") %>%
  mutate(
    p = ifelse(p == 0, 1e-322, p),
    log10P = -log10(p)
  ) %>%
  left_join(trait_category, by = "Trait")

category_order <- c(
  "Cardiovascular Diseases", "Diabetes Related", "Basic Traits",
  "Lipid","Blood Pressure", "Liver Function", "Renal Function",
  "Blood Cell", "Hormone Related", "Other"
)
score_df$Category <- factor(score_df$Category, levels = category_order)

score_df <- score_df %>%
  arrange(Category, desc(log10P)) %>%
  mutate(Trait = factor(Trait, levels = unique(Trait)))

brbg_colors <- c(
  "#003c30", "#543005",
  "#01665e", "#8c510a",
  "#35978f", "#bf812d",
  "#80cdc1", "#dfc27d",
  "#c7eae5", "#f6e8c3"
)

pdf("~/mgcv_maps_simulation/figures/manhatten_quantile.pdf",width=11,height=5.5)
ggplot(score_df, aes(x = Trait, y = log10P, color = Category)) +
  geom_segment(aes(xend = Trait, y = 0, yend = log10P), linewidth = 1.2) +
  geom_segment(aes(xend = Trait, y = 0, yend = log10P), linewidth = 0.2, color = "black") +
  geom_point(size = 2.8) +
  geom_point(size = 0.8, color = "black") +
  scale_color_manual(values = brbg_colors)+
  labs(
    x = "Trait",
    y = expression(-log[10](italic(P)))
  ) +
  geom_text(
    data = score_df[score_df$log10P > -log10(0.05/50), ],
    aes(x = Trait, y = log10P + 10, label = "*"),
    color = "black",
    size = 5,
    vjust = 0.5
  )+
  guides(color = guide_legend(title = NULL)) + 
  guides(color = guide_legend(title = NULL, nrow = 1)) + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),               
    legend.key.size = unit(0.5, "lines"),              
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)
  )
dev.off()
