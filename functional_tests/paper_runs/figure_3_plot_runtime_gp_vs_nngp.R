# Read the data from the file
library(ggplot2)
runtime_file <- "runtime_data_gp_vs_nngp.txt"
runtime_data <- read.table(runtime_file, col.names = c("seed",
                                                       "n",
                                                       "Method",
                                                       "runtime"))

# Calculate slope
slope_res <- coef(lm(log10(runtime) ~ log10(n) * Method, data = runtime_data))

g <- ggplot(runtime_data, aes(x = n, y = runtime, color = Method)) +
  geom_smooth(se = F, method = "lm") +
  geom_jitter(width = 0.005, height = 0.005, size = 1, alpha = 0.4) +
  scale_color_manual(values = c("GP" = "#4876B0", "nnGP" = "#B52A2B")) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  ylab("Wall Clock Time (s)") +
  xlab("Sample size") +
  theme_bw() +
  annotate(geom = "text", x = 5000, y = 200,
           label = sprintf("Slope = %.2f", slope_res[2]),
           color = "#4876B0") +
  annotate(geom = "text", x = 5000, y = 20,
           label = sprintf("Slope = %.2f", slope_res[2] + slope_res[4]),
           color = "#B52A2B")


plot(g)
