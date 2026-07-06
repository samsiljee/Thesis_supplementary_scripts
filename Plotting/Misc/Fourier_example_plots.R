# Code to make two example plots to illustrate Fourier transformation

data_point_resolution <- 1000
wave_1_amplitude <- 1
wave_1_frequency <- 1
wave_2_amplitude <- 0.7
wave_2_frequency <- 4

x <- seq(0,8*pi,length.out=data_point_resolution)

y_1 <- sin(wave_1_frequency * x) * wave_1_amplitude
y_2 <- sin(wave_2_frequency * x) * wave_2_amplitude
y_combined <- y_1 + y_2

data.frame(
  x = rep(x,3),
  y = c(y_1, y_2 + 1.8, y_combined + 5.5),
  colour = rep(c("red", "blue", "black"), each = data_point_resolution)
) %>%
  ggplot(aes(x = x, y = y, col = colour)) +
  geom_line(size = 1, show.legend = FALSE) +
  labs(
    x = "Time",
    y = element_blank()
  ) +
  scale_color_manual(values = c("black", "darkgreen", "darkblue")) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank())

# Save the plot
ggsave(paste0(output_dir, "Fourier waves.png"), width = 4, height = 5.5)

# Plot Fourier transform
data.frame(x = c(wave_1_frequency, wave_2_frequency), y = c(wave_1_amplitude, wave_2_amplitude), colour = c("darkgreen", "darkblue"))  %>%
  ggplot(aes(x = x, y = y, fill = colour)) +
  geom_col(width = 0.1, show.legend = FALSE) +
  scale_fill_manual(values = c("darkgreen", "darkblue")) +
  scale_colour_manual(values = c("darkgreen", "darkblue")) +
  theme_bw() +
  xlim(c(0, 5)) +
  labs(
    x = "Frequency",
    y = "Strength"
  ) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    )

# Save the plot
ggsave(paste0(output_dir, "Fourier transform.png"), width = 4, height = 5.5)
