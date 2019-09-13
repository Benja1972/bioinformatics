import pybedtools as pb
import sys,os
import math
import numpy as np

from Enh_Mut_Manip import bigWig2bed_pot

# Project settings
from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')


# H3K27ac signals
tlx_27ac = join(DATADIR,'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_H3K27ac_mm9.bw')
rag_27ac = join(DATADIR,'tracks/ChiP-seq_tracks/mm9_bigWig/RAG_H3K27ac_mm9.bw')



# TSS

tss = pb.BedTool(join(DATADIR,'tracks/annot_tracks/references/mm9/mm9.refGene.bed'))    # [1]

pad = 3e3 # 3k padding 

tss_tlx_27ac_3K = bigWig2bed_pot(tlx_27ac,tss,'mm9',pad=pad)

tss_tlx_27ac_3K.to_csv(join(DATADIR,'tracks/TSS_TLX3_H3K27ac_RP_3K.csv'))



#~ rag_27ac_tlx = bigWig2bed_pot(rag_27ac,enh,'mm9',pad=pad)

#~ enh_27ac_tlx.to_csv(join(DATADIR,'tracks/Enhancers_RAG_ChromHMM_TLX3_H3K27ac.csv'))
#~ rag_27ac_tlx.to_csv(join(DATADIR,'tracks/Enhancers_RAG_ChromHMM_RAG_H3K27ac.csv'))


# PolII 
#~ tlx_polII = join(DATADIR,'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_POLII_mm9.bw')

#~ enh_PolII_tlx = bigWig2bed_pot(tlx_polII,enh,'mm9',pad=pad)

#~ enh_PolII_tlx.to_csv(join(DATADIR,'tracks/Enhancers_ChromHMM_TLX3_PolII.csv'))






#enh =  enh.filter(lambda x: len(x) > 5900).saveas('temp.bed')

