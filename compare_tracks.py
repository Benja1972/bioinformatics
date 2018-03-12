#!/usr/bin/env python

import numpy as np
import yaml
from os.path import join 
import pybedtools as pb
import gffutils
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import matplotlib.pyplot as plt
import metaseq

#import subprocess as sb
#~ import svg_stack as ss
#~ import deeptools.bigwigCompare as bwComp
#~ import deeptools.multiBigwigSummary as bwCorr
#~ import deeptools.plotCorrelation as pltCorr

# make database from GTF
#~ db = gffutils.create_db("mouse.mm9.NCBIM37.67.ORIGINAL.gtf", "mouse.mm9.NCBIM37.67.ORIGINAL.gtf.db")

ref_dir = '../../references/mm9'
data_dir = 'tracks'
rel_path = "/home/sergio/media"

db = gffutils.FeatureDB(join(ref_dir,'mouse.mm9.NCBIM37.67.ORIGINAL.gtf.db'))


file_t = open(join(data_dir,"TLX3_TLX3_list.yaml"), "r")
file_r = open(join(data_dir,"RAG_TLX3_list.yaml"), "r")

tlx_lst = yaml.load(file_t)
rag_lst = yaml.load(file_r)



def tss_generator():
    """
    Generator function to yield TSS of each annotated transcript
    """
    for transcript in db.features_of_type('transcript'):
        yield TSS(asinterval(transcript), upstream=1, downstream=0)

# A BedTool made out of a generator, and saved to file.
#~ tsses = pb.BedTool(tss_generator()).saveas(join(ref_dir,'mm9_tsses.gtf'))

#Let's use pybedtools to add 1kb to either side of the TSS. This uses the BEDTools slop 
#routine; see the docs for that program for how to make changes to up/downstream distances.
#~ tsses_1kb = tsses.slop(b=1000, genome='mm9', output=join(ref_dir,'mm9_tsses-1kb.gtf'))
#~ tsses_2kb = tsses.slop(b=2000, genome='mm9', output=join(ref_dir,'mm9_tsses-2kb.gtf'))
#~ tsses_3kb = tsses.slop(b=3000, genome='mm9', output=join(ref_dir,'mm9_tsses-3kb.gtf'))
#~ tsses_5kb = tsses.slop(b=5000, genome='mm9', output=join(ref_dir,'mm9_tsses-5kb.gtf'))

tsses_1kb = pb.BedTool(join(ref_dir,'mm9_tsses-1kb.gtf'))


# ---------------------------------------------------------
#~ Metaseq #
# ---------------------------------------------------------

#~ metaseq works with the concepts of signal and windows. In this example, the signal is ChIP data, and the windows are TSS +/- 1kb.

#~ The first step is to create "genomic signal" objects out of the data. Since our example files are BAM files, we specify the kind='bam', but if you have your own data in a different format (bigWig, bigBed, BED, GFF, GTF, VCF) then specify that format instead (see metaseq.genomic_signal()).

#~ tlx_bdg = rel_path+tlx_lst["tracks"][0]
#~ rag_bdg = rel_path+rag_lst["tracks"][0]

tlx_bw = rel_path+tlx_lst["tracks"][2]
rag_bw = rel_path+rag_lst["tracks"][2]

print tlx_bw, rag_bw



tlx_signal = metaseq.genomic_signal(tlx_bw,'bigwig')
rag_signal = metaseq.genomic_signal(rag_bw,'bigwig')



