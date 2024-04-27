library(ggplot2)
runtime_file <- "runtime_data_var_nn_2.txt"
runtime_data <- NULL
runtime_data <- read.table(runtime_file, col.names = c("n",
                                                       "Method",
                                                       "nn_size",
                                                       "runtime"))
runtime_data$nn_size <- as.factor(runtime_data$nn_size)

# Define colors for each nn_size category
colors <- c("#B52A2B", "#609E3A", "#4876B0")

# Calculate slope for each nn_size
slopes <- by(runtime_data, runtime_data$nn_size, function(subset) {
  coef(lm(log10(runtime) ~ log10(n), data = subset))[2]
})

# Plot
g <- ggplot(runtime_data, aes(x = n, y = runtime, color = nn_size)) +
  geom_smooth(se = F, method = "lm") +
  geom_point(size = 1, alpha = 0.4) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  ylab("Wall clock time (s)") +
  xlab("Sample size") +
  theme_bw() +
  scale_color_manual(values = colors) # Apply manual colors

# Add slope annotations
for (i in seq_along(slopes)) {
  slope <- slopes[[i]]
  nn_size <- names(slopes)[i]
  g <- g +
    annotate(geom = "text", x = 8500, y = 10^i,
             label = sprintf("Slope (%s) = %.2f", nn_size, slope),
             color = colors[i])
}

plot(g)
