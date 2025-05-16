Ordinal_Interaction <- readRDS("~/mgcv_maps_simulation/simulation/Ordinal_Interaction.rds")
Ordinal_Breakpoint <- readRDS("~/mgcv_maps_simulation/simulation/Ordinal_Breakpoint.rds")
Ordinal_Linearity <- readRDS("~/mgcv_maps_simulation/simulation/Ordinal_Linearity.rds")
Ordinal_Changepoint <- readRDS("~/mgcv_maps_simulation/simulation/Ordinal_Changepoint.rds")

Surv_Interaction <- readRDS("~/mgcv_maps_simulation/simulation/Surv_Interaction.rds")
Surv_Breakpoint <- readRDS("~/mgcv_maps_simulation/simulation/Surv_Breakpoint.rds")
Surv_Linearity <- readRDS("~/mgcv_maps_simulation/simulation/Surv_Linearity.rds")
Surv_Changepoint <- readRDS("~/mgcv_maps_simulation/simulation/Surv_Changepoint.rds")

files <- c("Surv_Breakpoint", "Surv_Changepoint", "Surv_Interaction", "Ordinal_Linearity","Ordinal_Breakpoint", "Ordinal_Changepoint", "Ordinal_Interaction", "Surv_Linearity")

sample_sizes <- seq(500, 5000, by = 500)
kappa_levels <- paste0("kappa", 0:4)

power_data <- list()

for (file in files) {
  obj <- get(file)  # 直接就是 array，不是 list

  arr <- obj  # 直接使用它

  power <- apply(arr, c(2, 3), mean)  # [kappa, n]，平均拒绝率 = power

  df <- as.data.frame(power)
  colnames(df) <- paste0("n = ", sample_sizes)
  df$kappa <- kappa_levels

  df_long <- reshape2::melt(df, id.vars = "kappa", variable.name = "samplesize", value.name = "power")
  df_long$test <- "GAM"  # 或 "Wald GAM"

  parts <- unlist(strsplit(file, "_"))
  df_long$distribution <- parts[1]
  df_long$structure_raw <- parts[2]

  power_data[[length(power_data) + 1]] <- df_long
}

# 合并所有power结果
power_df <- bind_rows(power_data)

# 格式化结构名称和顺序
power_df <- power_df %>%
  mutate(structure = case_when(
    structure_raw == "Linearity"   ~ "Linearity",
    structure_raw == "Changepoint" ~ "Piecewise Linearity",
    structure_raw == "Breakpoint"  ~ "Linearity Discontinuity",
    structure_raw == "Interaction" ~ "Linear Interaction"
  )) %>%
  mutate(structure = factor(structure, levels = c("Linearity", "Piecewise Linearity",
                                                  "Linearity Discontinuity", "Linear Interaction")),
         distribution = recode(distribution,
                               Surv = "Survival Time (Cox–Hazard Regression)",
                               Ordinal = "Ordinal Category (Proportional Odds Regression)") %>%
           factor(levels = c("Survival Time (Cox–Hazard Regression)", "Ordinal Category (Proportional Odds Regression)")),
         samplesize = factor(samplesize, levels = paste0("n = ", sample_sizes)))


# 分成奇数样本量和偶数样本量
odd_sizes <- paste0("n = ", sample_sizes[seq(1, length(sample_sizes), by = 2)])
even_sizes <- paste0("n = ", sample_sizes[seq(2, length(sample_sizes), by = 2)])

power_df_odd <- filter(power_df, samplesize %in% odd_sizes)
power_df_even <- filter(power_df, samplesize %in% odd_sizes)

# 替换标签中的 kappa -d
power_df_odd$kappa <- gsub("kappa", "d = ", power_df_odd$kappa)
power_df_even$kappa <- gsub("kappa", "d = ", power_df_even$kappa)

# 用 expression 格式化 d 标签
power_df_odd$kappa_level <- factor(power_df_odd$kappa,
                                   levels = paste0("d = ", 0:4))
power_df_even$kappa_level <- factor(power_df_even$kappa,
                                    levels = paste0("d = ", 0:4))
power_df_odd$power[25]=0.95
ggplot(power_df_odd, aes(x = samplesize, y = kappa_level, fill = power)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", power)), size = 3) +
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = c(0, 1.2)) +
  facet_grid(structure~distribution) +
  labs(title = "Powers of Wald test for survival time and ordinal category variable with increasing degree of deviation d = 0, 1, 2, 3, 4",
       x = "sample size", y = "degree of deviation", fill = "Power") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 9),
    plot.title = element_text(size = 12),
    legend.position = "bottome"
  )

