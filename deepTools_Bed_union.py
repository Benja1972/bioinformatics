#!/usr/bin/env python


#~ import deeptools.bigwigCompare as bwComp
import deeptools.multiBigwigSummary as bwCorr
import deeptools.plotCorrelation as pltCorr
import numpy as np
import os
import yaml
#import subprocess as sb
import pybedtools as pb
import pandas as pd
import matplotlib.pyplot as plt
import svg_stack as ss

# plot style
#~ plt.style.use('ggplot')
#plt.style.use('fivethirtyeight')
#plt.style.use('seaborn-paper')

#~ params = {'font.size': 18,
          #~ 'lines.linewidth': 3
         #~ }
#~ plt.rcParams.update(params)          # Plot parameters


samples ={
    "TAP_H3K4me1",
    "TAP_H3K4me2",
    "TAP_H3K4me3",
    "TAP_H3K9ac",
    "TAP_H3K9me3",
    "TAP_H3K27ac",
    "TAP_H3K27me3",
    "TAP_H3k36me3",
    "TAP_POLII",
    "TAP_TLX3"
}

# directory name
np = '05peak/'
bp = '06broad_peak/'
#~ TLX3_H3K4me1_repl1_peaks.broadPeak

print '\\begin{center}'
print '\\begin{tabular}{ | l | l | l | p{5cm} |}'
print '\hline'
print 'Sample name & Num. peaks repl1 & Num. peaks repl2 & Num. peaks merged \\\\ \hline  \hline'
print '\\multicolumn{4}{|c|}{Narrow peaks}  \\\\ \hline'

for sample_name in samples:

	# files for replicates
	repl1 = np + sample_name + '_repl1_FE.bw'
	repl2 = np + sample_name + '_repl2_FE.bw'
	
	# Optional, for corellation on BED files of peaks, two files merge in one
	bed_repl1 = np + sample_name +'_repl1_summits.bed'
	bed_repl2 = np + sample_name +'_repl2_summits.bed'
	bed_merged = np + sample_name + '_merged.bed'
	bd_rp1 = pb.BedTool(bed_repl1)
	bd_rp2 = pb.BedTool(bed_repl2)
	bd_mrg = bd_rp1.cat(bd_rp2)
	bd_mrg.slop(l=1000,r=1000, genome = 'mm10').merge().saveas(bed_merged)
	msg = "\\verb|{}| & {} & {} & {} \\\\ \hline".format(sample_name, str(bd_rp1.count()),str(bd_rp2.count()), str(bd_mrg.count()))
	print msg



print '\hline \\multicolumn{4}{|c|}{Broad peaks}  \\\\ \hline'


for sample_name in samples:

	# files for replicates
	repl1 = bp + sample_name + '_repl1_FE.bw'
	repl2 = bp + sample_name + '_repl2_FE.bw'
	
	# Optional, for corellation on BED files of peaks, two files merge in one
	bed_repl1 = bp + sample_name +'_repl1_peaks.broadPeak' #'_repl1_summits.bed'
	bed_repl2 = bp + sample_name +'_repl2_peaks.broadPeak'  #'_repl2_summits.bed'
	bed_merged = bp + sample_name + '_merged.bed'
	bd_rp1 = pb.BedTool(bed_repl1)
	bd_rp2 = pb.BedTool(bed_repl2)
	bd_mrg = bd_rp1.cat(bd_rp2)
	bd_mrg.slop(l=1000,r=1000, genome = 'mm10').merge().saveas(bed_merged)
	msg = "\\verb|{}| & {} & {} & {} \\\\ \hline".format(sample_name, str(bd_rp1.count()),str(bd_rp2.count()), str(bd_mrg.count()))
	print msg

print '\end{tabular}'
print '\end{center}'

