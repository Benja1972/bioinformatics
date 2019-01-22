# Working with pdf
for f in *.eps; do epstopdf "$f"; done
pdftk *.pdf output TLX3-TAP-FE_states7_21.pdf

----------------------------------------------------------------------------------
# Enrichment by ChromHMM
ChromHMM.sh -Xmx4096M OverlapEnrichment TAP_6_FE_segments.bed enrichment TAP_enrichment > log 2>&1

---------------------------------------------------------------------------------
# GFF to BED
gff2bed < fimo.gff > RUNX_motif.bed

--------------------------------------------------------------------------
# bedtools slop will increase the size of each feature in a feature file by a user-defined number of bases

slopBed -i TLX3_TLX3_peaks.bed -g mm9 -b 100 > TLX3_TLX3_peaks_100.bed

--------------------------------------------------------------------------
# get DAN seqs from bed coordinates
bedtools getfasta -fi /home/sergio/media/NAS4/PFlab/TLX3_project/ChiP-Seq/references/mm9/chromFa/mm9.fa -bed TLX3_TLX3_peaks_100.bed -fo TLX3_TLX3_peaks_100.fa

 -----------------------------------------------------------------------------
# Extract specific state (ex. E7) from .bed

awk '($4=="E7")' file_segments.bed > file_E7_segments.bed

# for two states with OR
awk '($4=="E13" || $4=="E14")' file_segments.bed > file_E13E14_segments.bed


___________________________________________________________________
## LifOver
liftOver RP_mm10.bed mm10ToMm9.over.chain.gz RP_mm9.bed unmapped


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
# Remove 5th column in bed

cat file_in.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > file_out.bed