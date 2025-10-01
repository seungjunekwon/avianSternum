# Load libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Read the data with gene symbols (1st column) and normalized counts from DESeq2 (2nd column~)
dat <- read.csv("dat.csv") 

# Reshape the data from wide format to long format
data_long <- dat %>%
  pivot_longer(cols = starts_with("e"),
               names_to = "condition",
               values_to = "Normalized_Expression") %>%
  mutate(
    # Create a new variable to categorize the embryonic stage
    estage = case_when(
      grepl("e24", condition) ~ "St. 24",
      grepl("e36d", condition) ~ "St. 36_D",
      grepl("e36v", condition) ~ "St. 36_V"
    ),
    # Replicate the gene symbols for each condition
    symbol = rep(dat$symbol, each = 9)
  )

# Convert the 'symbol' column to a factor to ensure ggplot uses it correctly
data_long$symbol <- factor(data_long$symbol, levels = unique(dat$symbol))

# Define color palette for the conditions
st_colors <- c("St. 24" = "#0033CC", "St. 36_D" = "#00CC66", "St. 36_V" = "#FF69B4")

# Create the box plot with jitter for each gene
a <- ggplot(data_long, aes(x = estage, y = Normalized_Expression, color = estage)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.5) +
  geom_jitter(position = position_jitter(0.2), size = 2, alpha = 0.7) +
  scale_y_continuous(
    limits = c(0, 150),
    breaks = c(0, 50, 100, 150)
  ) +
  scale_color_manual(values = st_colors) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(color = "black", size = 9),
    legend.title = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold.italic", hjust = 0.5, margin = margin(b = 30)),
    panel.spacing = unit(2, "lines")
  ) +
  labs(y = "Normalized expression") +
  guides(color = guide_legend(title = "Stage"))

a
ggsave("boxPlot.png", plot = a, width = 8, height = 5, dpi = 300, units = "cm")
