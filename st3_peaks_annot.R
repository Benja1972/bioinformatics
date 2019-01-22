#source("https://bioconductor.org/biocLite.R")
#biocLite(c("ChIPpeakAnno", "ChIPseeker"))

#biocLite(c("TxDb.Mmusculus.UCSC.mm9.knownGene","EnsDb.Mmusculus.v75", "org.Mm.eg.db"))

#library(ChIPpeakAnno)
library(ChIPseeker)

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library("org.Mm.eg.db")
library(EnsDb.Mmusculus.v75)


smpl = "st3_peaks.bed"

track <- readPeakFile(smpl)

txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

track_anno <-annotatePeak(track, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Mm.eg.db")

genes_df = track_anno@anno@elementMetadata


setEPS()
postscript(paste(smpl,"_annotPie.eps",sep=""), horizontal = FALSE)
plotAnnoPie(track_anno)
dev.off()

postscript(paste(smpl,"_upset.eps",sep=""), width=8,height=4, horizontal = FALSE)
upsetplot(track_anno, vennpie=TRUE)
dev.off()
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

gene <- seq2gene(track, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)


pathway <- enrichPathway(gene, organism = "mouse")

#postscript(paste(smpl,"_pathway.eps", sep=""), paper="special", width = 650, height = 408, horizontal = FALSE)
postscript(paste(smpl,"_pathway.eps", sep=""),width=7,height=4, horizontal = FALSE)
dotplot(pathway)
dev.off()


# genes out 
genes2 = unique(track_anno@anno@elementMetadata@listData$ENSEMBL)
write.csv(genes2,paste(smpl,"_genes.csv", sep=""))

## ======= ChiPeakAnno for transcrips ======

#source("https://bioconductor.org/biocLite.R")
#biocLite(c("ChIPpeakAnno", "ChIPseeker", "toGRanges"))

library(ChIPpeakAnno)


trackA <- toGRanges(smpl , format="BED",  header=FALSE)

annoData <-toGRanges(EnsDb.Mmusculus.v75, feature="transcript")

trackA_anno <- annotatePeakInBatch(trackA,
                                     AnnotationData=annoData)


genesT = trackA_anno@elementMetadata@listData$feature
# write.csv(near_genesT,paste(smpl,"_genesT.csv", sep=""))




## ========== Genes expression analysis ======= 
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("DESeq2", "pheatmap","dplyr"))


library(DESeq2)

counts <- read.table('TLX3vsRAG_featureCounts.txt', header=TRUE, row.names=1)

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

nt <- normTransform(dds) # defaults to log2(x+1)
nt = assay(nt)

nt_cut = subset(nt, rownames(nt) %in% genesT)

select <- order(rowMeans(nt_cut),decreasing=TRUE) #[1:40]
nt_cut = nt_cut[select,]
df <- as.data.frame(colData(dds))


postscript(paste(smpl,"_geneExpr.eps", sep=""),width=7,height=14, horizontal = FALSE)
pheatmap(nt_cut, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()

gene_names	=	read.table("UCSC_mm9_transcripID_to_geneSymbol.sort.txt", row.names=1)
gene_namesC = subset(gene_names, rownames(gene_names) %in%  rownames(nt_cut))


library(dplyr)


gene_namesC<- add_rownames(gene_namesC, "GENES")
nt_cutC<- add_rownames(as.data.frame(nt_cut), "GENES")

gene_namesCC = gene_namesC[match(nt_cutC$GENES,gene_namesC$GENES),]

#nt_cutCC = nt_cutC[match(gene_namesC$GENES,nt_cutC$GENES),]

write.csv(gene_namesCC,paste(smpl,"_gene_names.csv", sep=""))





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


