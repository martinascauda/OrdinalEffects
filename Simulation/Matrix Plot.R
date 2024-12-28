library(patchwork)
library(ggplot2)
library(ggrain)
load("./Simulation/DataInt17_For_Plots.RData")
data$Level = factor(data$Level, levels = sort(unique(data$Level)), ordered = TRUE)
data$Outcome = factor(data$Outcome, levels = sort(unique(data$Outcome)), ordered = TRUE)
data$OCE <- as.numeric(data$OCE)


# Create individual plots
plot_list <- lapply(unique(data$Outcome), function(outcome_value) {
  ggplot(data[data$Outcome == outcome_value, ], 
         aes(x = Method, y = OCE, fill = Level)) +
    # Delete the scatterplot to enhance readability of the plot
    geom_rain(alpha = 0.5, rain.side = 'f', # cov = "Level",
              boxplot.args = list(outlier.shape = NA, alpha = 0.8, linewidth=0.3),
              violin.args = list(alpha = 0.5, color = NA),
              point.args = list(size = 0.001, color="white"), #0.1
              point.args.pos = list(position = position_jitterdodge(
                jitter.width = 0.0001,
                jitter.height = 0,
                dodge.width = 0.5,
                seed = 42
              )),
              boxplot.args.pos = list(width = 0.15,
                                      position = ggpp::position_dodgenudge(width = 0.5,
                                                                           x = c(rep(-0.045, 6), rep(0.045, 6)))),
              violin.args.pos = list(width = 0.9,
                                     position = position_nudge(
                                       x = c(rep(rep(-0.3, 256 * 2), 6), rep(rep(0.3, 256 * 2), 6)))
              )) +
    scale_y_continuous(
      breaks = seq(round(min(data$OCE), digits=1), round(max(data$OCE), digits=1), by = 0.2) # Adjust break interval
    ) +
    stat_summary(fun = mean, geom = "line", aes(group = Level, color = Level),
                 position = position_dodge(width = 0.5), linewidth = 0.3, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", aes(group = Level, color = Level),
                 position = position_dodge(width = 0.5), size = 0.5, pch = 23, alpha = 1, color = "black") +
    scale_fill_manual(name = "Outcome Levels",
                      labels = 1:6, values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                                               "#66A61E", "#E6AB02", "#A6761D", "#666666")) +
    scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                                  "#66A61E", "#E6AB02", "#A6761D", "#666666")) +
    theme_classic() +
    labs(x = "", y = "OCE", 
         title = paste("Outcome Variable:", outcome_value)) +
    theme(
      legend.position = "none",  # Remove legend for individual plots
      plot.title = element_text(hjust = 0.5)  # Center plot titles
    )
})

# Create a shared legend as a separate plot
legend_plot <- ggplot(data, aes(x = Method, y = OCE, fill = Level)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(name = "Outcome\nLevels",
                    labels = 1:6, values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                                             "#66A61E", "#E6AB02", "#A6761D", "#666666")) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 10)
  )

# Extract the legend
legend = cowplot::get_plot_component(legend_plot, 'guide-box-top', return_all = TRUE)
legend <- cowplot::ggdraw(legend)

# Combine plots with the legend inside the grid
ncol <- 4  # Number of columns
nrow <- ceiling(length(plot_list) / ncol)  # Calculate rows based on plots and columns

# Convert legend into a patchwork plot
legend_as_plot <- wrap_elements(grid::grid.draw(legend)) 

# Add legend after the last plot in the grid
combined_plot <- wrap_plots(c(plot_list, list(legend_as_plot)), ncol = ncol)

ggsave(file="Simulation/Matrix17.pdf", width = 11.55, height = 10.69, dpi = 300)
