# Load libraries
library(ggplot2)
library(dplyr)
library(ggtext)

# Read the data
dat <- read.csv('go.csv') # 3 columns: 'GO BP hierarchy', 'foldEnrichment','P'; PANTHER hierarchy terms were used

# Compute -log10(P)
dat$negLog10P <- -log10(dat$P)

# Highlight terms
highlight_terms <- c("Extracellular matrix organization", 
                     "Positive regulation of cell population proliferation")

# HTML styling of GO terms
dat$GO.BP.label <- ifelse(dat$GO.BP.hierarchy %in% highlight_terms,
                          paste0("<b style='color:black;'>", dat$GO.BP.hierarchy, "</b>"),
                          paste0("<span style='color:gray40;'>", dat$GO.BP.hierarchy, "</span>"))

# Maintain display order
dat$GO.BP.label <- factor(dat$GO.BP.label,
                          levels = dat$GO.BP.label[order(dat$foldEnrichment, decreasing = FALSE)])

# Plot
a <- ggplot(dat, aes(x = foldEnrichment, y = GO.BP.label)) +
  geom_point(aes(size = negLog10P), shape = 21, fill = "#FF69B4", color = "black", stroke = 0.4) +
  theme_minimal(base_family = "Arial") +
  theme(
    # Axis styling with ggtext
    axis.text.y.right = element_markdown(
      hjust = 0,
      size = 10,
      margin = margin(l = 5)
    ),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(),
    
    # Hide left y-axis
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.title.y.left = element_blank(),
    
    # Hide grid, show full box
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    axis.line = element_blank(),
    
    # Legend placement
    legend.position = "left"
  ) +
  scale_size_continuous(
    name = expression(-log[10] ~ p),
    range = c(3, 8),  # 30% bigger dots
    guide = guide_legend(override.aes = list(
      shape = 21, fill = NA, color = "black", stroke = 0.4))
  ) +
  scale_x_continuous(name = "Fold Enrichment", limits=c(2,18)) +
  scale_y_discrete(position = "right") +
  coord_cartesian(clip = "off")

a

ggsave("go.png", plot = a, width = 14, height = 7, dpi = 300, units = "cm")
