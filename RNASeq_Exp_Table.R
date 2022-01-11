# Import libraries for RNA-Seq analysis and plots
library("DESeq2")
library("ggplot2")

# Read sample files within current folder, use grep to filter and return matched files as vector
sampleFiles <- grep("count", list.files(), value = TRUE)
sampleFiles
 
# Design the sample table dataframe, read.table with pheotypes of large sample size, do not hard coding
sampleNames <- c("BG07","BG08","BG09","BG10","BG11","BG12")
sampleBatch <- c("Batch1","Batch2","Batch3","Batch1","Batch2","Batch3")
sampleCondition <- c("WT_B","WT_B","WT_B","KO_B","KO_B","KO_B")

sampleTable <- data.frame(sampleName = sampleNames, 
                          fileName = sampleFiles, 
                          condition = sampleCondition, 
                          batch = sampleBatch)

# Assign the reference level (control group), otherwise use contrast for complex design
sampleTable$condition <- factor(sampleTable$condition, levels = c("WT_B","KO_B"))
sampleTable

# Combine the ht-seq count data and sample information to generate DESeq object
ddsHTseqCounts <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~ condition)
ddsHTseqCounts

# Pre-filtering the dataset to remove low/no expression gene
ddsHTseq <- ddsHTseqCounts[rowSums(counts(ddsHTseqCounts)) >= 10, ]
ddsHTseq

# Run the DESeq pipeline, estimating size factors, dispersions, fitting model and Wald test
dds <- DESeq(ddsHTseq)
dds

# Extracts a results table with log2 fold changes, p values and adjusted p values
# cooksCutoff is disabled to get unfiltered results, which contains outlier genes
res <- results(dds,cooksCutoff = FALSE)
res

# Order results by padj value (most significant to least)
resOrdered <- res[order(res$padj),]
summary(resOrdered)

# Merge the DESeq2 object results with normalized counts from each sample
# Use row.names (gene names) as the master key for merge 

resdata <- merge(as.data.frame(resOrdered), as.data.frame(counts(dds,normalized=T)), by='row.names', sort=F)

# Change the first cell from 'Row.names' to 'gene'

names(resdata)[1] <- 'gene'
head(resdata)

# Save results and normalized counts to files

write.csv(resdata, file="Culture_Normalized_GBM_B.csv")
