# Creates interactive volcano plot in Rstudio

library(ggplot2)
library(plotly)
library(readr)

# Load your data
df <- read_tsv("DE_results_annotated.tsv")

# Check column names (important!)
print(colnames(df))

# Define regulation categories
df$regulation <- "Not significant"
df$regulation[df$padj < 0.05 & df$log2FoldChange > 1]  <- "Upregulated"
df$regulation[df$padj < 0.05 & df$log2FoldChange < -1] <- "Downregulated"

# Plot with Annotation in hover text
p <- ggplot(df, aes(x = log2FoldChange, 
                    y = -log10(padj),
                    color = regulation,
                    text = paste(
                      "ID:", IDs,
                      "<br>Annotation:", Annotation,
                      "<br>log2FC:", round(log2FoldChange, 2),
                      "<br>padj:", signif(padj, 3),
                      "<br>baseMean:", round(baseMean, 2)
                    ))) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano plot with annotations",
       x = "log2(Fold Change)",
       y = "-log10(adjusted p-value)",
       color = "Regulation")

# Convert to interactive plot
interactive_plot <- ggplotly(p, tooltip = "text")
interactive_plot
