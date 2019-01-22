source("https://bioconductor.org/biocLite.R")
biocLite(c("DESeq", "pheatmap","dplyr"))


genesT = read.table('state6_tlx_not_rag_genesT.csv', header=TRUE, row.names=1, sep = ',')
genesT = as.matrix(genesT)

counts <- read.table('TLX3vsRAG_featureCounts.txt', header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
counts <- counts[ ,6:ncol(counts)]

# Convert to matrix
counts <- as.matrix(counts)

# Assign condition (first three are controls, second three contain the expansion)
(condition <- factor(c(rep("rag", 3), rep("tlx", 3))))

library(DESeq2)
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(counts), condition))

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)

# Run the DESeq pipeline
dds <- DESeq(dds)

# Regularized log transformation for clustering/heatmaps, etc
# rld <- rlogTransformation(dds)

# Get differential expression results
# res <- results(dds)

## Order by adjusted p-value
# res <- res[order(res$padj), ]
## Merge with normalized count data
# resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=TRUE)), by="row.names", sort=FALSE)
# names(resdata)[1] <- "Gene"

library("pheatmap")



# ts = counts(dds,normalized=TRUE)
# ts_cut = subset(ts, rownames(ts) %in% genesT)
# select <- order(rowMeans(ts_cut),decreasing=TRUE)[1:40]

nt <- normTransform(dds) # defaults to log2(x+1)
nt = assay(nt)

nt_cut = subset(nt, rownames(nt) %in% genesT)
select <- order(rowMeans(nt_cut),decreasing=TRUE)[1:40]
nt_cut = nt_cut[select,]
df <- as.data.frame(colData(dds))


pheatmap(nt_cut, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)





gene_names	=	read.table("UCSC_mm9_transcripID_to_geneSymbol.sort.txt", row.names=1)
gene_namesC = subset(gene_names, rownames(gene_names) %in%  rownames(nt_cut))


library(dplyr)


gene_namesC<- add_rownames(gene_namesC, "GENES")
nt_cutC<- add_rownames(as.data.frame(nt_cut), "GENES")

gene_namesCC = gene_namesC[match(nt_cutC$GENES,gene_namesC$GENES),]

#nt_cutCC = nt_cutC[match(gene_namesC$GENES,nt_cutC$GENES),]

write.csv(gene_namesCC,"genes40.csv")

# 
# 
# nt_cutt = subset(nt_cut, rownames(nt_cut) %in% rownames(gene_namesC))
# 
# #nt_cutt = order(rownames(gene_namesC))
# 
# nt_cutt = as.data.frame.data.frame(nt_cutt)
# gene_namesC = as.data.frame(gene_namesC)
# 
# 
# 
# library(dplyr)
# 
# 
# nt_cutt <- add_rownames(nt_cutt, "GENES")
# gene_namesC <- add_rownames(gene_namesC, "GENES")
# 
# 
# colnames(gene_namesC)[colnames(gene_namesC)=="V2"] <- "Genes_names"
# 
# all_dt = inner_join(nt_cutt,gene_namesC)


