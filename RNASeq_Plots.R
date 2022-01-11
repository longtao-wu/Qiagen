# The code for plots of RNA-Seq data, adjust axis limit according to expression table

# PCA Method 1: Apply a 'regularized log' transformation

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
plotPCA(rld)

p <- plotPCA(rld)
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
print(p)

# PCA Method 2: Variance stabilization and normalization

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

q <- plotPCA(vsd)
q <- q + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
print(q)


# Dispersion Plot

postscript("Dispersion_plot.eps", paper = "letter", colormodel = "rgb", onefile= FALSE, horizontal = FALSE, family = "Times")
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# MA Plots of different flavors

postscript("MA_plot_01_v2.eps", paper = "letter", colormodel = "rgb", onefile= FALSE, horizontal = FALSE, family = "Times")

plot(res$baseMean, res$log2FoldChange, pch=16, cex=0.5, log="x", ylim=c(-10,13), main="MA Plot")

  with(subset(res, padj<0.01), points(baseMean, log2FoldChange, col="yellow", pch=16, cex=0.5))
  with(subset(res, padj<0.01 & log2FoldChange>0), points(baseMean, log2FoldChange, col="red", pch=16, cex=0.5))
  with(subset(res, padj<0.01 & log2FoldChange<0), points(baseMean, log2FoldChange, col="green", pch=16, cex=0.5))

dev.off()

postscript("MA_plot_02_v2.eps", paper = "letter", colormodel = "rgb", onefile= FALSE, horizontal = FALSE, family = "Times")

plot(res$baseMean, res$log2FoldChange, pch=16, cex=0.5, log="x", ylim=c(-10,13), main="MA Plot")

  with(subset(res, padj<0.01), points(baseMean, log2FoldChange, col="yellow", pch=16, cex=0.5))
  with(subset(res, padj<0.01 & log2FoldChange>=1), points(baseMean, log2FoldChange, col="red", pch=16, cex=0.5))
  with(subset(res, padj<0.01 & log2FoldChange<=-1), points(baseMean, log2FoldChange, col="green", pch=16, cex=0.5))

dev.off()

postscript("MA_plot_03_v2.eps", paper = "letter", colormodel = "rgb", onefile= FALSE, horizontal = FALSE, family = "Times")

plot(res$baseMean, res$log2FoldChange, pch=16, cex=0.5, log="x", ylim=c(-10,13), main="MA Plot")

  with(subset(res, padj<0.01 & log2FoldChange>=1), points(baseMean, log2FoldChange, col="red", pch=16, cex=0.5))
  with(subset(res, padj<0.01 & log2FoldChange<=-1), points(baseMean, log2FoldChange, col="green", pch=16, cex=0.5))

dev.off()

# Volcano plot with genes pass minimal expression level, low expression dots mess up plots
res_mean_name <- rownames(subset(res, baseMean >=20 ))
res_mean <- res[res_mean_name,]

plot(res_mean$log2FoldChange, -log10(res_mean$padj), pch=16, main="Volcano Plot", ylim=c(0,25), xlim=c(-5,10))
  with(subset(res_mean, padj<0.01 ), points(log2FoldChange, -log10(padj), pch=16, col="yellow"))
  with(subset(res_mean, padj<0.01 & log2FoldChange>=1), points(log2FoldChange, -log10(padj), pch=16, col="red"))
  with(subset(res_mean, padj<0.01 & log2FoldChange<=-1), points(log2FoldChange, -log10(padj), pch=16, col="green"))

abline(h=2, col="orange", pch=32)
abline(v=c(-1,1), col="black", pch=32)