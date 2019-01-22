#source("https://bioconductor.org/biocLite.R")
#biocLite(c("ChIPpeakAnno", "ChIPseeker"))

#biocLite(c("TxDb.Mmusculus.UCSC.mm9.knownGene","EnsDb.Mmusculus.v75", "org.Mm.eg.db"))

#library(ChIPpeakAnno)
library(ChIPseeker)

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library("org.Mm.eg.db")
library(EnsDb.Mmusculus.v75)

txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

tlx_f="tracks/TLX3_TLX3_peaks_100.bed" 
tlx_tr <- readPeakFile(tlx_f)

tlx_anno <-annotatePeak(tlx_tr, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

df_anno <-as.data.frame(tlx_anno)

write.csv(df_anno,paste(substr(tlx_f,start=1, stop=26),"_annot.csv", sep=""))

#genes_df = track_anno@anno@elementMetadata

# genes out 
# tlx_gene = unique(tlx_anno@anno@elementMetadata@listData$ENSEMBL)

#write.csv(genes2,paste(smpl,"_genes.csv", sep=""))


## ======= ChiPeakAnno for transcrips ======

#source("https://bioconductor.org/biocLite.R")
#biocLite(c("ChIPpeakAnno", "ChIPseeker", "toGRanges"))

# library(ChIPpeakAnno)
# 
# annoData <-toGRanges(EnsDb.Mmusculus.v75, feature="transcript")
# tlx_tra <- toGRanges(tlx_f , format="BED",  header=FALSE)
# 
# tlx_annA <- annotatePeakInBatch(tlx_tra, AnnotationData=annoData)
# 
# 
# tlx_geneT = tlx_annA@elementMetadata@listData$feature
# 
# 
# write.csv(tlx_geneT,paste(tlx_f,"_genesT.csv", sep=""))



## ====== Annotation figures 


# setEPS()
# postscript(paste(tlx_f,"_annotPie.eps",sep=""), horizontal = FALSE)
plotAnnoPie(tlx_anno)
# dev.off()

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
# library(ReactomePA)
# 
# tlx_geneR <- seq2gene(tlx_tr, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
# 
# 
# tlx_pathway <- enrichPathway(tlx_geneR, organism = "mouse")
# 
# #postscript(paste(smpl,"_pathway.eps", sep=""), paper="special", width = 650, height = 408, horizontal = FALSE)
# postscript(paste(tlx_f,"_pathway.eps", sep=""),width=12,height=4, horizontal = FALSE)
# dotplot(tlx_pathway)
# dev.off()


