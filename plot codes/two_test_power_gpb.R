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
kappa_levels <- paste0("kappa", 0:4)

power_data <- list()

for (file in files) {
  obj <- get(file)
  for (test_type in c("B1", "B2")) {
    arr <- obj[[test_type]]  

    power <- apply(arr, c(2, 3), mean)  

    df <- as.data.frame(power)
    colnames(df) <- paste0("n = ", sample_sizes)
    df$kappa <- kappa_levels

    df_long <- reshape2::melt(df, id.vars = "kappa", variable.name = "samplesize", value.name = "power")
    df_long$test <- ifelse(test_type == "B1", "Wald test", "Score test")

    parts <- unlist(strsplit(file, "_"))
    df_long$distribution <- parts[1]
    df_long$structure_raw <- parts[2]

    power_data[[length(power_data) + 1]] <- df_long
  }
}
power_df <- bind_rows(power_data)
power_df <- power_df %>%
  mutate(structure = case_when(
    structure_raw == "Linearity"   ~ "Linearity",
    structure_raw == "Changepoint" ~ "Piecewise Linearity",
    structure_raw == "Breakpoint"  ~ "Linearity discontounity",
    structure_raw == "Interaction" ~ "Linear Interaction"
  )) %>%
  mutate(structure = factor(structure, levels = c("Linearity", "Piecewise Linearity",
                                                  "Linearity discontounity", "Linear Interaction")),
         distribution = factor(distribution, levels = c("Gaussian", "Binary", "Poisson")),
         samplesize = factor(samplesize, levels = paste0("n = ", sample_sizes)))

odd_sizes <- paste0("n = ", sample_sizes[seq(1, length(sample_sizes), by = 2)])
even_sizes <- paste0("n = ", sample_sizes[seq(2, length(sample_sizes), by = 2)])

power_df_odd <- filter(power_df, samplesize %in% odd_sizes)
power_df_even <- filter(power_df, samplesize %in% odd_sizes)

power_df_odd$kappa <- gsub("kappa", "d = ", power_df_odd$kappa)
power_df_even$kappa <- gsub("kappa", "d = ", power_df_even$kappa)

power_df_odd$kappa_level <- factor(power_df_odd$kappa,
                                   levels = paste0("d = ", 0:4))
power_df_even$kappa_level <- factor(power_df_even$kappa,
                                    levels = paste0("d = ", 0:4))

power_df_odd=power_df_odd[which(power_df_odd$test=="Score test"),]
power_df_even=power_df_even[which(power_df_even$test=="Wald test"),]

plot3 = ggplot(power_df_odd, aes(x = samplesize, y = kappa_level, fill = power)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", power)), size = 3) +
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = c(0, 1.2)) +
  facet_grid(distribution~structure) +
  labs(title = "C. Powers of score test with increasing degree of deviation d = 0, 1, 2, 3, 4",
       x = "sample size", y = "degree of deviation", fill = "Power") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 9),
    plot.title = element_text(size = 12),
    legend.position = "right"
  )

plot4 = ggplot(power_df_even, aes(x = samplesize, y = kappa_level, fill = power)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", power)), size = 3) +
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = c(0, 1.2)) +
  facet_grid(distribution~structure) +
  labs(title = "D. Powers of Wald test with increasing degree of deviation d = 0, 1, 2, 3, 4",
       x = "sample size", y = "degree of deviation", fill = "Power") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 9),
    plot.title = element_text(size = 12),
    legend.position = "right"
  )

egg::ggarrange(plot1,plot2,plot3,plot4)
