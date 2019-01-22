#source("https://bioconductor.org/biocLite.R")
#biocLite(c("ChIPpeakAnno", "ChIPseeker"))

#biocLite(c("TxDb.Mmusculus.UCSC.mm9.knownGene","EnsDb.Mmusculus.v75", "org.Mm.eg.db"))

#library(ChIPpeakAnno)
library(ChIPseeker)

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library("org.Mm.eg.db")
library(EnsDb.Mmusculus.v75)

txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

tlx_f="tracks/state17_PolII_badRAG/TLX_st16only.bed"
rag_f="tracks/state17_PolII_badRAG/RAG_st16only.bed"
tap_f="tracks/state17_PolII_badRAG/TAP_st16only.bed"


tlx_tr <- readPeakFile(tlx_f)
rag_tr <- readPeakFile(rag_f)
tap_tr <- readPeakFile(tap_f)




tlx_anno <-annotatePeak(tlx_tr, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Mm.eg.db")
rag_anno <-annotatePeak(rag_tr, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Mm.eg.db")
tap_anno <-annotatePeak(tap_tr, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Mm.eg.db")

#genes_df = track_anno@anno@elementMetadata





# genes out 
tlx_gene = unique(tlx_anno@anno@elementMetadata@listData$ENSEMBL)
rag_gene = unique(rag_anno@anno@elementMetadata@listData$ENSEMBL)
tap_gene = unique(tap_anno@anno@elementMetadata@listData$ENSEMBL)


#write.csv(genes2,paste(smpl,"_genes.csv", sep=""))


## ======= ChiPeakAnno for transcrips ======

#source("https://bioconductor.org/biocLite.R")
#biocLite(c("ChIPpeakAnno", "ChIPseeker", "toGRanges"))

library(ChIPpeakAnno)

annoData <-toGRanges(EnsDb.Mmusculus.v75, feature="transcript")

tlx_tra <- toGRanges(tlx_f , format="BED",  header=FALSE)
rag_tra <- toGRanges(rag_f , format="BED",  header=FALSE)
tap_tra <- toGRanges(tap_f , format="BED",  header=FALSE)



tlx_annA <- annotatePeakInBatch(tlx_tra, AnnotationData=annoData)
rag_annA <- annotatePeakInBatch(rag_tra, AnnotationData=annoData)
tap_annA <- annotatePeakInBatch(tap_tra, AnnotationData=annoData)


tlx_geneT = tlx_annA@elementMetadata@listData$feature
rag_geneT = rag_annA@elementMetadata@listData$feature
tap_geneT = tap_annA@elementMetadata@listData$feature


# write.csv(near_genesT,paste(smpl,"_genesT.csv", sep=""))



## ====== Annotation figures 


setEPS()
postscript(paste(tlx_f,"_annotPie.eps",sep=""), horizontal = FALSE)
plotAnnoPie(tlx_anno)
dev.off()

setEPS()
postscript(paste(rag_f,"_annotPie.eps",sep=""), horizontal = FALSE)
plotAnnoPie(rag_anno)
dev.off()

setEPS()
postscript(paste(tap_f,"_annotPie.eps",sep=""), horizontal = FALSE)
plotAnnoPie(tap_anno)
dev.off()

# postscript(paste(smpl,"_upset.eps",sep=""), width=8,height=4, horizontal = FALSE)
# upsetplot(track_anno, vennpie=TRUE)
# dev.off()
#upsetplot(track02_Anno vennpie=TRUE)

# ======= Functional enrichment analysis ======= 
# 
# Once we have obtained the annotated nearest genes, we can perform functional enrichment 
# analysis to identify predominant biological themes among these genes by incorporating biological
# knowledge provided by biological ontologies. For instance, Gene Ontology (GO)7 annotates genes 
# to biological processes, molecular functions, and cellular components in a directed acyclic 
# graph structure, Kyoto Encyclopedia of Genes and Genomes (KEGG)8 annotates genes to pathways, 
# Disease Ontology (DO)9 annotates genes with human disease association, and Reactome10 annotates 
# gene to pathways and reactions.
# 
# ChIPseeker also provides a function, seq2gene, for linking genomc regions to genes in 
# a many-to-many mapping. It consider host gene (exon/intron), promoter region and flanking 
# gene from intergenic region that may under control via cis-regulation. This function is 
# designed to link both coding and non-coding genomic regions to coding genes and facilitate functional analysis.
# 
# Enrichment analysis is a widely used approach to identify biological themes. I have developed 
# several Bioconductor packages for investigating whether the number of selected genes 
# associated with a particular biological term is larger than expected, including DOSE2 
# for Disease Ontology, ReactomePA for reactome pathway, clusterProfiler4 for Gene Ontology 
# and KEGG enrichment analysis.
library(ReactomePA)

tlx_geneR <- seq2gene(tlx_tr, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
rag_geneR <- seq2gene(rag_tr, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
tap_geneR <- seq2gene(tap_tr, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)


tlx_pathway <- enrichPathway(tlx_geneR, organism = "mouse")
rag_pathway <- enrichPathway(rag_geneR, organism = "mouse")
tap_pathway <- enrichPathway(tap_geneR, organism = "mouse")

#postscript(paste(smpl,"_pathway.eps", sep=""), paper="special", width = 650, height = 408, horizontal = FALSE)
postscript(paste(tlx_f,"_pathway.eps", sep=""),width=7,height=4, horizontal = FALSE)
dotplot(tlx_pathway)
dev.off()

postscript(paste(rag_f,"_pathway.eps", sep=""),width=7,height=4, horizontal = FALSE)
dotplot(rag_pathway)
dev.off()

postscript(paste(tap_f,"_pathway.eps", sep=""),width=7,height=4, horizontal = FALSE)
dotplot(tap_pathway)
dev.off()





## ========== Genes expression analysis ======= 
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("DESeq2", "pheatmap","dplyr"))


library(DESeq2)

counts <- read.table('TLX3vsRAG_featureCounts.txt', header=TRUE, row.names=1)
gene_names = read.table("UCSC_mm9_transcripID_to_geneSymbol.sort.txt", row.names=1)

# Remove first five columns (chr, start, end, strand, length)
counts <- counts[ ,6:ncol(counts)]

# Convert to matrix
counts <- as.matrix(counts)

# Assign condition (first three are controls, second three contain the expansion)
(condition <- factor(c(rep("rag", 3), rep("tlx", 3))))

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(counts), condition))


dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)

# Run the DESeq pipeline
dds <- DESeq(dds)

library("pheatmap")
library(dplyr)

nt <- normTransform(dds) # defaults to log2(x+1)
nt = assay(nt)

tlx_cut = subset(nt, rownames(nt) %in% tlx_geneT)
rag_cut = subset(nt, rownames(nt) %in% rag_geneT)
tap_cut = subset(nt, rownames(nt) %in% tap_geneT)

tlx_sel <- order(rowMeans(tlx_cut),decreasing=TRUE) #[1:40]
rag_sel <- order(rowMeans(rag_cut),decreasing=TRUE)[1:120]
tap_sel <- order(rowMeans(tap_cut),decreasing=TRUE) #[1:40]

tlx_cut = tlx_cut[tlx_sel,]
rag_cut = rag_cut[rag_sel,]
tap_cut = tap_cut[tap_sel,]


df <- as.data.frame(colData(dds))


postscript(paste(tlx_f,"_geneExpr.eps", sep=""),width=7,height=14, horizontal = FALSE)
pheatmap(tlx_cut, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off()

postscript(paste(rag_f,"_geneExpr.eps", sep=""),width=7,height=14, horizontal = FALSE)
pheatmap(rag_cut, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off()

postscript(paste(tap_f,"_geneExpr.eps", sep=""),width=7,height=14, horizontal = FALSE)
pheatmap(tap_cut, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off()


tlx_gene_namesC = subset(gene_names, rownames(gene_names) %in%  rownames(tlx_cut))
rag_gene_namesC = subset(gene_names, rownames(gene_names) %in%  rownames(rag_cut))
tap_gene_namesC = subset(gene_names, rownames(gene_names) %in%  rownames(tap_cut))





tlx_gene_namesC<- add_rownames(tlx_gene_namesC, "GENES")
rag_gene_namesC<- add_rownames(rag_gene_namesC, "GENES")
tap_gene_namesC<- add_rownames(tap_gene_namesC, "GENES")

tlx_cutC<- add_rownames(as.data.frame(tlx_cut), "GENES")
rag_cutC<- add_rownames(as.data.frame(rag_cut), "GENES")
tap_cutC<- add_rownames(as.data.frame(tap_cut), "GENES")

tlx_gene_namesCC = tlx_gene_namesC[match(tlx_cutC$GENES,tlx_gene_namesC$GENES),]
rag_gene_namesCC = rag_gene_namesC[match(rag_cutC$GENES,rag_gene_namesC$GENES),]
tap_gene_namesCC = tap_gene_namesC[match(tap_cutC$GENES,tap_gene_namesC$GENES),]



write.csv(tlx_gene_namesCC,paste(tlx_f,"_gene_names.csv", sep=""))
write.csv(rag_gene_namesCC,paste(rag_f,"_gene_names.csv", sep=""))
write.csv(tap_gene_namesCC,paste(tap_f,"_gene_names.csv", sep=""))





#horiz=TRUE,onefile=FALSE,width=8.5,height=11,paper=letter)
# ==========================================================================
# Find peaks with bi-directional promoters
# 
# Bidirectional promoters are the DNA regions located between TSS of two adjacent genes that are transcribed 
# on opposite directions and often co-regulated by this shared promoter region5. Here is an example to find 
# peaks near bi-directional promoters.

# bdp <- peaksNearBDP(tlx_gr, annoData, maxgap=5000)

# c(bdp$percentPeaksWithBDP, 
#   bdp$n.peaks, 
#   bdp$n.peaksWithBDP)





# covplot(tlx_e6)
# 
# promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=500, by = "transcript")
# tagMatrix <- getTagMatrix(tlx_e6, windows=promoter)
# 
# tagHeatmap(tagMatrix, xlim=c(-2000, 500), color="red")


