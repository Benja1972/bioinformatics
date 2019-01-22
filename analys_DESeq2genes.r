# source("https://bioconductor.org/biocLite.R")
# biocLite(c("edgeR"))

library(DESeq2)
#library(yaml)


tb <- "TLX3vsRAGvsTAP_featCounts_genes.txt"
counts <- read.table(tb, header=TRUE, row.names=1)
counts <- counts[ ,6:ncol(counts)]
counts <- as.matrix(counts)

condition <- factor(c(rep("rag", 3), rep("tlx", 3), rep("tap", 3)))
coldata <- data.frame(row.names=colnames(counts), condition)

# Run the DESeq pipeline
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- DESeq(dds)


#plotDispEsts(dds, main="Dispersion plot")

## === Get differential expression results ===
# res <- results(dds)

# === TLX3 vs RAG ===
res_TLX_RAG <- results(dds, contrast=c("condition","tlx","rag"))
#table(res$padj<0.05)
## Order by adjusted p-value
res_TLX_RAG <- res_TLX_RAG[order(res_TLX_RAG$padj), ]
## Merge with normalized count data
res_TLX_RAGout <- merge(as.data.frame(res_TLX_RAG), as.data.frame(counts(dds, normalized=TRUE)[,1:6]), by="row.names", sort=FALSE)
names(res_TLX_RAGout)[1] <- "Gene"

## Write results
write.table(res_TLX_RAGout, file="TLX3vsRAG-results_genes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# === TAP vs RAG ===
res_TAP_RAG <- results(dds, contrast=c("condition","tap","rag"))
#table(res$padj<0.05)
## Order by adjusted p-value
res_TAP_RAG <- res_TAP_RAG[order(res_TAP_RAG$padj), ]
## Merge with normalized count data
res_TAP_RAGout <- merge(as.data.frame(res_TAP_RAG), as.data.frame(counts(dds, normalized=TRUE)[,c(1,2,3,7,8,9)]), by="row.names", sort=FALSE)
names(res_TAP_RAGout)[1] <- "Gene"

## Write results
write.table(res_TAP_RAGout, file="TAPvsRAG-results_genes.txt", row.names = FALSE, quote = FALSE, sep = "\t")


# === TAP vs TLX3
res_TAP_TLX3 <- results(dds, contrast=c("condition","tap","tlx"))

## Order by adjusted p-value
res_TAP_TLX3 <- res_TAP_TLX3[order(res_TAP_TLX3$padj), ]
## Merge with normalized count data
res_TAP_TLX3out <- merge(as.data.frame(res_TAP_TLX3), as.data.frame(counts(dds, normalized=TRUE)[,4:9]), by="row.names", sort=FALSE)
names(res_TAP_TLX3out)[1] <- "Gene"

## Write results
write.table(res_TAP_TLX3out, file="TAPvsTLX3-results_genes.txt", row.names = FALSE, quote = FALSE, sep = "\t")


# cooks = assays(dds)[["cooks"]]
###################### ===================== ##################################################################

#resdata["ENSMUST00000103276", "RAGS.RAGZ"]

#res <- results(dds)
#hist(res$pvalue, breaks=50, col="grey")

# library(gplots)
# library(RColorBrewer)
# library("pheatmap")
# 
# # Regularized log transformation for clustering/heatmaps, etc
# rld <- rlogTransformation(dds)
# mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]
# # Sample distance heatmap
# sampleDists <- as.matrix(dist(t(assay(rld))))
# heatmap.2(as.matrix(sampleDists), key=F, trace="none",
#           col=colorpanel(100, "black", "white"),
#           ColSideColors=mycols[condition], RowSideColors=mycols[condition],
#           margin=c(10, 10), main="Sample Distance Matrix")
# 
# plotPCA(rld, intgroup=c("condition"))
# 
# 
# DESeq2::plotMA(dds,  ylim=c(-15,15), cex=1)


# volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
#   with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
#   with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
#   with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
#   with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
#   if (labelsig) {
#     require(calibrate)
#     #with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
#   }
#   legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
# }
# #png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
# volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-15, 15))

#dev.off()

# select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:55]
# nt <- normTransform(dds) # defaults to log2(x+1)
# df <- as.data.frame(colData(dds))
# 
# pheatmap(assay(nt)[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)

# ===== EdgeR =====
# library(edgeR)
# 
# y <- DGEList(counts=counts, group=condition)
