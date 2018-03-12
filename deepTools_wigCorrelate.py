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
plt.style.use('ggplot')
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
bw = '05peak/'
brd_bw = '06broad_peak/'

for sample_name in samples:
	# Sample name for replicates correlation
	#sample_name = 'TLX3_TLX3'


	# files for replicates
	repl1 = bw + sample_name + '_repl1_FE.bw'
	repl2 = bw + sample_name + '_repl2_FE.bw'
	
	# Optional, for corellation on BED files of peaks, two files merge in one
	bed_merged = brd_bw + sample_name + '_merged.bed'
	
	# labels for graphs
	lb1 = sample_name + '_r1'
	lb2 = sample_name + '_r2'

	# output file as compressed numpy array
	out_file = bw + sample_name + '_replCorr' + '.npz'
	out_fileB = bw + sample_name + '_replCorr' + '.npz'

	# output figure in SVG format
	out_figP  = bw + sample_name + '_replCorr_pear' + '.svg'
	out_figS  = bw + sample_name + '_replCorr_spear' + '.svg'
	out_fig  = bw + sample_name + '_replCorr' + '.svg'   

	print repl1, repl2, out_file

	# --binSize 50
	#args = "bins -b {} {}  -o {}".format(repl1, repl2, out_file).split()
	args = "BED-file -b {} {}  -o {} --binSize 500 --BED {}".format(repl1, repl2, out_fileB, bed_merged).split()
	bwCorr.main(args)
	#~ if not(os.path.exists(out_file)):  
		#~ args = "bins -b {} {}  -o {}".format(repl1, repl2, out_file).split()
		#~ #args = "BED-file -b {} {}  -o {} --BED {}".format(repl1, repl2, out_file, bed_merged).split()
		#~ bwCorr.main(args)

	# Pearson correlation
	argsP = "-in {} --whatToPlot scatterplot --labels {} {} --corMethod {} -o {}".format(out_fileB, lb1, lb2, 'pearson', out_figP).split()
	pltCorr.main(argsP)

	# Spearman correlation
	argsS = "-in {} --whatToPlot scatterplot --labels {} {} --corMethod {} -o {}".format(out_fileB, lb1, lb2, 'spearman', out_figS).split()
	pltCorr.main(argsS)
    
    #~ # Combine two fugures in one file with svg_stack utility
	#~ doc = ss.Document()

	#~ # layout = ss.VBoxLayout()
	#~ layout = ss.HBoxLayout()
	#~ layout.addSVG(out_figP,alignment=ss.AlignTop|ss.AlignHCenter)
	#~ layout.addSVG(out_figS,alignment=ss.AlignCenter)

	#~ doc.setLayout(layout)

	#~ doc.save(out_fig)

# Combine all fugures in one file with svg_stack utility
docAll = ss.Document()
layoutA = ss.VBoxLayout()


for sample_name in samples:
	#~ out_fig  = bw + sample_name + '_replCorr' + '.svg'
	out_figS  = bw + sample_name + '_replCorr_spear' + '.svg'
	layoutA.addSVG(out_figS,alignment=ss.AlignCenter)

docAll.setLayout(layoutA)

docAll.save(bw+sample_name[:4]+'_all_replCorr' + '.svg')





