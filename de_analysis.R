# Load libraries
library(tximport)
#library(readr)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# === 1. Read metadata ===
samples <- read.csv("samples.csv")
files <- file.path(paste0("quant_", samples$sample), "quant.sf")
names(files) <- samples$sample

# === 2. Transcript → gene mapping ===
tx2gene <- read.table("tx2gene.tsv", sep="\t", header=FALSE)
colnames(tx2gene) <- c("TXNAME","GENEID")

# === 3. Import counts ===
txi <- tximport(files, type="salmon", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM",
		ignoreTxVersion=TRUE)

# === 4. Build DESeq2 object ===
dds <- DESeqDataSetFromTximport(txi,
            colData=samples,
            design=~condition)

# === 5. Filter low-count genes ===
dds <- dds[rowSums(counts(dds)) > 10, ]

# === 6. Run DESeq2 ===
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","treatment","control"))
res <- lfcShrink(dds, contrast=c("condition","treatment","control"),
                 type="ashr")   # shrink log2FCs for stability
summary(res)

# === 7. Save results ===
resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered), file="DE_results.csv")

# === 8. QC plots ===
# PCA
vsd <- vst(dds)
pca <- plotPCA(vsd, intgroup="condition")
ggsave("PCA_plot.png", pca)

# Volcano
res_df <- as.data.frame(resOrdered)
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(alpha=0.5) +
  geom_vline(xintercept=c(-1,1), col="red", lty=2) +
  geom_hline(yintercept=-log10(0.05), col="red", lty=2) +
  theme_minimal() + ggtitle("Volcano plot")

# Top-variable genes heatmap
topgenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 50)
annotation_col <- data.frame(condition=samples$condition)
rownames(annotation_col) <- samples$sample  # ensure alignment
pheatmap(assay(vsd)[topgenes, ],
         annotation_col=annotation_col,
         show_rownames=FALSE,
         cluster_cols=TRUE)

