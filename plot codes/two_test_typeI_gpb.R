library(dplyr)
library(reshape2)
library(ggplot2)

Binary_Breakpoint <- readRDS("~/mgcv_maps_simulation/simulation/Binary_Breakpoint.rds")
Binary_Interaction <- readRDS("~/mgcv_maps_simulation/simulation/Binary_Interaction.rds")
Poisson_Breakpoint <- readRDS("~/mgcv_maps_simulation/simulation/Poisson_Breakpoint.rds")
Poisson_Interaction <- readRDS("~/mgcv_maps_simulation/simulation/Poisson_Interaction.rds")
Binary_Changepoint <- readRDS("~/mgcv_maps_simulation/simulation/Binary_Changepoint.rds")
Binary_Linearity <- readRDS("~/mgcv_maps_simulation/simulation/Binary_Linearity.rds")
Poisson_Linearity <- readRDS("~/mgcv_maps_simulation/simulation/Poisson_Linearity.rds")
Poisson_Changepoint <- readRDS("~/mgcv_maps_simulation/simulation/Poisson_Changepoint.rds")
Gaussian_Interaction <- readRDS("~/mgcv_maps_simulation/simulation/Gaussian_Interaction.rds")
Gaussian_Linearity <- readRDS("~/mgcv_maps_simulation/simulation/Gaussian_Linearity.rds")
Gaussian_Changepoint <- readRDS("~/mgcv_maps_simulation/simulation/Gaussian_Changepoint.rds")
Gaussian_Breakpoint <- readRDS("~/mgcv_maps_simulation/simulation/Gaussian_Breakpoint.rds")

files <- c("Binary_Breakpoint", "Binary_Changepoint", "Binary_Interaction", "Binary_Linearity",
           "Gaussian_Breakpoint", "Gaussian_Changepoint", "Gaussian_Interaction", "Gaussian_Linearity",
           "Poisson_Breakpoint", "Poisson_Changepoint", "Poisson_Interaction", "Poisson_Linearity")

sample_sizes <- seq(500, 5000, by = 500)

all_data <- list()

for (file in files) {
  obj <- get(file)
  for (test_type in c("P1", "P2")) {
    mat <- obj[[test_type]]
    df <- as.data.frame(mat)
    colnames(df) <- sample_sizes
    df$sim_id <- 1:nrow(df)
    df_long <- melt(df, id.vars = "sim_id", variable.name = "samplesize", value.name = "pvalue")
    df_long$samplesize <- paste0("n = ",as.numeric(as.character(df_long$samplesize)))
    df_long$test <- ifelse(test_type == "P1", "Wald test", "Score test")

    parts <- unlist(strsplit(file, "_"))
    df_long$distribution <- parts[1]
    df_long$structure_raw <- parts[2]

    all_data[[length(all_data) + 1]] <- df_long
  }
}

plot_data <- bind_rows(all_data)

plot_data <- plot_data %>%
  mutate(structure = case_when(
    structure_raw == "Linearity"   ~ "Linearity",
    structure_raw == "Changepoint" ~ "Piecewise Linearity",
    structure_raw == "Breakpoint"  ~ "Linearity discontounity"
  )) %>%
  mutate(structure = factor(structure, levels = c("Linearity", "Piecewise Linearity",
                                                  "Linearity discontounity", "Linear Interaction")),
         distribution = factor(distribution, levels = c("Binary", "Poisson", "Gaussian")))

plot_data_sorted <- plot_data %>%
  group_by(distribution, structure, samplesize, test) %>%
  arrange(pvalue) %>%
  mutate(p_sorted = pvalue,
         x_theory = (1:n()) / n())

plot_data_sorted$distribution=ordered(plot_data_sorted$distribution,levels=c("Gaussian","Binary","Poisson"))
plot_data_sorted$samplesize=ordered(plot_data_sorted$samplesize,levels=paste0("n = ",seq(500,5000,500)))
library(RColorBrewer)

greens_palette <- brewer.pal(9, "BrBG")[c(1:3,8:9)]

odd_sizes <- paste0("n = ", seq(500, 4500, by = 1000))  
even_sizes <- paste0("n = ", seq(1000, 5000, by = 1000))

plot_data_odd <- plot_data_sorted %>%
  filter(test == "Score test", samplesize %in% odd_sizes)

plot_data_even <- plot_data_sorted %>%
  filter(test == "Wald test", samplesize %in% odd_sizes)

plot1 <- ggplot(plot_data_odd, aes(x = x_theory, y = p_sorted, color = samplesize, group = samplesize)) +
  geom_abline(slope = 1, intercept = 0, color = "gray80", size = 3) +
  geom_line(size = 1.5) +
  scale_color_manual(values = greens_palette, name = "Sample size") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  facet_grid(distribution ~ structure) +
  labs(
    title = "A. QQ plot of p-values of score test when the null hypothesis test holds",
    x = "uniform quantile",
    y = "sorted p-values",
    color = "Test Method"
  ) +
  theme(strip.text = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.position = "bottom")

plot2 <- ggplot(plot_data_even, aes(x = x_theory, y = p_sorted, color = samplesize, group = samplesize)) +
  geom_abline(slope = 1, intercept = 0, color = "gray80", size = 3) +
  geom_line(size = 1.5) +
  scale_color_manual(values = greens_palette, name = "Sample size") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  facet_grid(distribution ~ structure) +
  labs(
    title = "B. QQ plot of p-values of Wald test when the null hypothesis test holds",
    x = "uniform quantile",
    y = "sorted p-values",
    color = "Test Method"
  ) +
  theme(strip.text = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.position = "bottom")
