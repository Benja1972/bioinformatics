#source("https://bioconductor.org/biocLite.R")
#biocLite(c("ChIPpeakAnno", "ChIPseeker"))

#biocLite(c("TxDb.Mmusculus.UCSC.mm9.knownGene","EnsDb.Mmusculus.v75", "org.Mm.eg.db"))

#library(ChIPpeakAnno)
library(ChIPseeker)

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library("org.Mm.eg.db")
library(EnsDb.Mmusculus.v75)

# = RAG
#smpl = "RAG_14_E6"
#smpl = "RAG_14_E13"

# = TLX3
#smpl = "TLX3_14_E6"
#smpl = "TLX3_14_E13"

# = complex
#smpl = "state6_tlx_not_rag"
smpl = "RAG_st3"

#file_nm = paste(smpl,"_sorted.bed", sep="")
file_nm = paste(smpl,".bed", sep="")

# =============== ChiPeakAnno
# track01 <- toGRanges(file_nm , format="BED",  header=FALSE)
# 
# annoData <-toGRanges(EnsDb.Mmusculus.v75, feature="gene")
# annoDataT <-toGRanges(EnsDb.Mmusculus.v75, feature="transcript")

# track01_anno <- annotatePeakInBatch(track01,
#                                      AnnotationData=annoData,
#                                      output="nearestBiDirectionalPromoters",
#                                      bindingRegion=c(-2000, 500))

# track01_anno <- annotatePeakInBatch(track01,
#                                     AnnotationData=annoData)
# 
# track01_annoT <- annotatePeakInBatch(track01,
#                                     AnnotationData=annoDataT)

# nearest genes
# near_genes = track01_anno@elementMetadata@listData$feature
# write.csv(near_genes,paste(smpl,"_genes.csv", sep=""))

# near_genesT = track01_annoT@elementMetadata@listData$feature
# write.csv(near_genesT,paste(smpl,"_genesT.csv", sep=""))

# genes
# genes01 = track01_anno@elementMetadata@listData$gene_name
# genes_sorted = sort(genes01)
# write.csv(genes_sorted,paste(smpl,"_genesBI.csv", sep=""))

# =============== ChIPseeker
track02 <- readPeakFile(file_nm)

txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

track02_anno <-annotatePeak(track02, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Mm.eg.db")

genes_df = track02_anno@anno@elementMetadata


#track02_anno <-annotatePeak(track02, tssRegion=c(-1000, 1000), TxDb=txdb)
#genes02 = track02_anno@anno@elementMetadata@listData$SYMBOL
#genes_sorted02 = sort(genes02)
#write.csv(genes_sorted02,paste(smpl,"_genesLong.csv", sep=""))


setEPS()
postscript(paste(smpl,"_annotPie.eps",sep=""), horizontal = FALSE)
plotAnnoPie(track02_anno)
dev.off()

postscript(paste(smpl,"_upset.eps",sep=""), width=8,height=4, horizontal = FALSE)
upsetplot(track02_anno, vennpie=TRUE)
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

gene03 <- seq2gene(track02, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)


pathway <- enrichPathway(gene03, organism = "mouse")

#postscript(paste(smpl,"_pathway.eps", sep=""), paper="special", width = 650, height = 408, horizontal = FALSE)
postscript(paste(smpl,"_pathway.eps", sep=""),width=7,height=4, horizontal = FALSE)
dotplot(pathway)
dev.off()


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


