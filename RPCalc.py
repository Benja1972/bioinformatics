import pybedtools as pb
import sys,os
import math
import numpy as np
import deeptools.getScorePerBigWigBin as gs






def bigWig2bed_pot(bw,bed,genome,pad=1e5,step =1, alpha=10):
    bed_df = bed.to_dataframe()
    bed_df['mid'] = (bed_df['end'] + bed_df['start'])/2
    bed_df['mid'] = bed_df['mid'].astype('int')
    bed_df['mid_end'] = bed_df['mid'] + 1

    col = ['chrom','mid','mid_end', 'name']

    tss = pb.BedTool.from_dataframe(bed_df[col])

    pad = int(pad)
    z = np.arange(-pad,pad+1,step)
    wt = 2.0*np.exp(-alpha*np.fabs(z)/1e5)/(1.0+np.exp(-alpha*np.fabs(z)/1e5))

    df = tss.slop(b=pad, genome=genome).to_dataframe()

    for i, row in df.iterrows():

        cnt = tss[i].start
        st = df.loc[i,'start']
        end = df.loc[i,'end']
        chrom = df.loc[i,'chrom']

        #~ cmd = '{0} {1} {2} {3} {4} {5}'.format('bigWigSummary',bw, chrom, st, end, end-st)

        #~ result = os.popen(cmd)
        #~ content = result.readlines()
        #~ temp = content[0].strip().replace('n/a','0')
        #~ temp = temp.split('\t')
        #~ values = np.array(temp,dtype=np.float64)
        #print(chrom, st, end, [bw])
        
        
        vl = gs.countFragmentsInRegions_worker(chrom, int(st), int(end), [bw], stepSize=step,binLength=step, save_data=False)
        vl = np.transpose(np.squeeze(vl[0]))
        vl  = np.hstack((np.zeros(st - cnt + pad),vl,np.zeros(cnt + pad - end +1 )))

        df.loc[i,'RP'] = np.dot( vl, wt)
    
    return df


# Project settings
from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')

### Main

# Enh TLX3
#~ enh = pb.BedTool(join(DATADIR,'tracks/Enhancers_ChromHMM.bed'))

# Enh RAG
enh = pb.BedTool(join(DATADIR,'tracks/Enhancers_RAG_ChromHMM.bed'))

pad = 3e3 # 3k padding 

# H3K27ac
tlx_27ac = join(DATADIR,'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_H3K27ac_mm9.bw')
rag_27ac = join(DATADIR,'tracks/ChiP-seq_tracks/mm9_bigWig/RAG_H3K27ac_mm9.bw')

enh_27ac_tlx = bigWig2bed_pot(tlx_27ac,enh,'mm9',pad=pad)
rag_27ac_tlx = bigWig2bed_pot(rag_27ac,enh,'mm9',pad=pad)

enh_27ac_tlx.to_csv(join(DATADIR,'tracks/Enhancers_RAG_ChromHMM_TLX3_H3K27ac.csv'))
rag_27ac_tlx.to_csv(join(DATADIR,'tracks/Enhancers_RAG_ChromHMM_RAG_H3K27ac.csv'))


# PolII 
#~ tlx_polII = join(DATADIR,'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_POLII_mm9.bw')

#~ enh_PolII_tlx = bigWig2bed_pot(tlx_polII,enh,'mm9',pad=pad)

#~ enh_PolII_tlx.to_csv(join(DATADIR,'tracks/Enhancers_ChromHMM_TLX3_PolII.csv'))






#enh =  enh.filter(lambda x: len(x) > 5900).saveas('temp.bed')

