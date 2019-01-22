
# Extract specific states (ex. E7) from .bed

awk '($4=="E13")' TLX3_14_segments.bed | sort-bed - > TLX3_st3.bed
awk '($4=="E6")' TLX3_14_segments.bed | sort-bed - > TLX3_st8.bed
awk '($4=="E5")' TLX3_14_segments.bed | sort-bed - > TLX3_st9.bed

awk '($4=="E13")' RAG_14_segments.bed | sort-bed - > RAG_st3.bed
awk '($4=="E6")' RAG_14_segments.bed | sort-bed - > RAG_st8.bed
awk '($4=="E5")' RAG_14_segments.bed | sort-bed - > RAG_st9.bed


# Combine beds
cat TLX3_st8.bed TLX3_st9.bed | sort-bed - | uniq > TLX3_st8-9.bed
cat RAG_st8.bed RAG_st9.bed | sort-bed - | uniq > RAG_st8-9.bed


__________________________________________
# Sort bed

sort -k1,1 -k2,2n in.bed > in.sorted.bed

sort-bed unsorted.TSS.bed > TSS.bed

# combined command
awk '($4=="E13")' TLX3_14_segments.bed | sort-bed - > TLX3_14_E13_sorted.bed
______________________________________________
# Intersect beds
Compute the overlap between your transcription factor peaks and the UCSC promoters using bedtools by typing the following command line:
bedtools (http://bedtools.readthedocs.org/en/latest/) is a set of tools to manipulate BED files. Here, we will use the bedtools intersect (or intersectBed) tool to compute the overlap between two BED files.

intersectBed -a XXX -b YYY -wa -wb -f 0.10 | uniq > ZZZ.inter

Replace XXX by the filename of your BED file
Replace YYY by the filename of the UCSC promoters/enhancers BED file
Replace ZZZ by the output filename (e.g. GATA3_overlap_promoters)

Parameters of bedtools intersect tool are available by typing intersectBed --help in a UNIX Terminal :
-a is followed by the 1st BED file
-b is followed by the 2nd BED file
-wa indicates that we want the regions of the 1st file which overlap with regions of the 2nd file
-f indicates that we want a minimum overlap of 10% (0.10) for the regions in the 1st file

#
closest-features TLX3_14_E6_sorted.bed genes.bed > TLX3_E6_genes.bed

---------------------------------------
# Collect all genes

$ awk '{print $8}' TLX3_E6_overlap_genes.inter | uniq > TLX3_E6_genes.txt


Looking at your conditions, it sounds like you always want the closest downstream element, regardless of what's upstream or how near or far away that upstream element is.

By default, the BEDOPS closest-features application will report the full characteristics of both the nearest leftmost ("upstream") and rightmost ("downstream") elements.

You can therefore feed the output from this application to awk (or another interpreted language) to just report the TSS element and whatever nearest gene is downstream from it.

As a for instance, first sort your input files:

$ sort-bed unsorted.TSS.bed > TSS.bed
$ sort-bed unsorted.genes.bed > genes.bed
Then run the two sorted BED files through closest-features, piping the results to awk with the correct field separator to pick the downstream gene:

$ closest-features TSS.bed genes.bed \
    | awk FS="|" '{ \
        tssElement = $1; \
        upstreamGene = $2; \
        downstreamGene = $3; \  
        print tssElement"|"downstreamGene; \
    }' - \
    > answer.bed
The file answer.bed will be a sorted BED file that contains results in the following format:

[ tss-1 ] | [ nearest-downstream-gene-to-tss-1 ]
[ tss-2 ] | [ nearest-downstream-gene-to-tss-2 ]
...
[ tss-n ] | [ nearest-downstream-gene-to-tss-n ]
If you have a different condition for picking a nearest element, you could set up those rules in the awk block above. (Feel free to modify your question if you had something else in mind.)

To test features on the basis of strand, you could modify the awk block as follows:

$ closest-features TSS.bed genes.bed \
    | awk FS="|" '{ \
        tssRegion = $1; \
        upstreamGene = $2; \
        downstreamGene = $3; \  
        split(tssRegion, tssElements, "\t"); \
        split(upstreamGene, upstreamGeneElements, "\t"); \
        split(downstreamGene, downstreamGeneElements, "\t"); \
        tssStart = tssElements[1]; \
        tssStop = tssElements[2]; \
        upstreamGeneStart = upstreamGeneElements[1]; \
        upstreamGeneStop = upstreamGeneElements[2]; \
        upstreamGeneStrand = upstreamGeneElements[4]; \
        downstreamGeneStart = downstreamGeneElements[1]; \
        downstreamGeneStop = downstreamGeneElements[2]; \
        downstreamGeneStrand = downstreamGeneElements[4]; \
        trueUpstreamGeneDistance = 0; \
        trueDownstreamGeneDistance = 0; \
        if (upstreamGeneStrand == "+") { \
            trueUpstreamGeneDistance = tssStart - upstreamGeneStart; \
        } \
        else { \
            trueUpstreamGeneDistance = tssStart - upstreamGeneStop; \
        } \
        if (downstreamGeneStrand == "+") { \
            trueDownstreamGeneDistance = downstreamGeneStart - tssStop; \
        } \
        else { \
            trueDownstreamGeneDistance = downstreamGeneStop - tssStop; \
        } \
        if (trueUpstreamGeneDistance > trueDownstreamGeneDistance) { \
            print tssRegion"|"downstreamGene; \
        } \
        else { \
            print tssRegion"|"upstreamGene; \
        } \
    }' - \
    > answer.bed