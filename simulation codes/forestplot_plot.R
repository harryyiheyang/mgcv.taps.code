G1 <- readRDS("~/mgcv_maps_simulation/RDD/AgeBreak6070.rds")
G2 <- readRDS("~/mgcv_maps_simulation/RDD/AgeChange6070.rds")
G1$jump.q=FDRestimation::p.fdr(G1$jump.p,adjust.method = "BY")$fdrs
G2$slopechange.q=FDRestimation::p.fdr(G2$slopechange.p,adjust.method = "BY")$fdrs
G11=G1[which(G1$Score.p>0.05&G1$jump.q<0.1),]
G22=G2[which(G2$Score.p>0.05&G2$slopechange.q<0.1),]
G11=G11%>%dplyr::select(Outcome=outcome,`Jump Estimate`=jump.est,`Jump SE`=`jump.se`,`Jump P`=jump.p,`Jump Q`=jump.q,`Score Test P`=Score.p)
G22=G22[-11,]

#library(forestplot)
#library(dplyr)

G11 <- G1 %>%
  filter(Score.p > 0.05 & jump.q < 0.1) %>%
  transmute(
    Outcome = outcome,
    `Jump Estimate` = jump.est,
    `Jump SE` = jump.se,
    `Jump P` = jump.p,
    `Jump Q` = jump.q,
    `Score Test P` = Score.p
  ) %>%
  arrange(Outcome,`Jump Q`)

# ---------- 构建绘图数据框 ----------
G11 %>%
  mutate(
    mean = `Jump Estimate`,
    lower = `Jump Estimate` - 1.96 * `Jump SE`,
    upper = `Jump Estimate` + 1.96 * `Jump SE`,
    `Score P` = formatC(`Score Test P`, format = "e", digits = 2),
    `Jump P` = formatC(`Jump P`, format = "e", digits = 2),
    `Jump Q` = formatC(`Jump Q`, format = "e", digits = 2),
    `Jump Estimate` = sprintf("%.3f", `Jump Estimate`)
  ) |>
  forestplot(
  labeltext = c("Outcome", "Score P", "Jump P", "Jump Q"),
  is_summary = FALSE,
  xlab = "Jump Estimate (95% CI)",
  box_size = 0.25,
  lineheight = unit(16, "pt"),
  xticks = c(-0.1,-0.05, 0, 0.05),
  ci_column = TRUE
)  |>
  fp_add_lines() |>
  fp_set_style(box = "black",
               line = "black",
               summary = "#60c5ba")  |>
  fp_add_header(Outcome="Outcome",
                `Score P` = "Score Test P",
                `Jump P` = "Jump P",
                `Jump Q` = "Jump Q") |>
  fp_set_zebra_style("#EFEFEF")
###############################################################################
G11 <- G22 %>%
  filter(Score.p > 0.05 & slopechange.q < 0.1) %>%
  transmute(
    Outcome = outcome,
    `Jump Estimate` = slopechange.est,
    `Jump SE` = slopechange.se,
    `Jump P` = slopechange.p,
    `Jump Q` = slopechange.q,
    `Score Test P` = Score.p
  ) %>%
  arrange(Outcome,`Jump Q`)

# ---------- 构建绘图数据框 ----------
G11 %>%
  mutate(
    mean = `Jump Estimate`,
    lower = `Jump Estimate` - 1.96 * `Jump SE`,
    upper = `Jump Estimate` + 1.96 * `Jump SE`,
    `Score P` = formatC(`Score Test P`, format = "e", digits = 2),
    `Jump P` = formatC(`Jump P`, format = "e", digits = 2),
    `Jump Q` = formatC(`Jump Q`, format = "e", digits = 2),
    `Jump Estimate` = sprintf("%.3f", `Jump Estimate`)
  ) |>
  forestplot(
    labeltext = c("Outcome", "Score P", "Jump P", "Jump Q"),
    is_summary = FALSE,
    xlab = "Slope Change Estimate (95% CI)",
    box_size = 0.25,
    lineheight = unit(16, "pt"),
    xticks = c(-0.1,-0.05, 0, 0.05),
    ci_column = TRUE
  )  |>
  fp_add_lines() |>
  fp_set_style(box = "black",
               line = "black",
               summary = "#60c5ba")  |>
  fp_add_header(Outcome="Outcome",
                `Score P` = "Score Test P",
                `Jump P` = "Slope P",
                `Jump Q` = "Slope Q") |>
  fp_set_zebra_style("#EFEFEF")
