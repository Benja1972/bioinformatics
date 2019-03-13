import pybedtools as pb
import sys,os
import math
import numpy as np

# Project settings
from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')

#al = 1e4
#alpha = -math.log(1.0/3.0)*1e5/al
alpha=10
bw = join(DATADIR,'tracks/ChiP-seq_tracks/mm10_bigWig/TAP_H3K27me3_treat_pileup.bw')

enh = pb.BedTool(join(DATADIR,'tracks/Enhancers_ChromHMM.bed'))
bd = pb.BedTool(join(DATADIR,'tracks/UP_RP_enhancers.bed'))

def bigWig2bed_pot(bw,bed,genome,pad=1e5,alpha=10):
    bed_df = bed.to_dataframe()
    bed_df['mid'] = (bed_df['end'] + bed_df['start'])/2
    bed_df['mid'] = bed_df['mid'].astype('int')
    bed_df['mid_end'] = bed_df['mid'] + 1

    col = ['chrom','mid','mid_end', 'name']

    tss = pb.BedTool.from_dataframe(bed_df[col])

    pad = int(pad)
    z = np.arange(-pad,pad+1)
    wt = 2.0*np.exp(-alpha*np.fabs(z)/1e5)/(1.0+np.exp(-alpha*np.fabs(z)/1e5))

    df = tss.slop(b=pad, genome=genome).to_dataframe()

    for i, row in df.iterrows():

        cnt = tss[i].start
        st = df.loc[i,'start']
        end = df.loc[i,'end']
        chrom = df.loc[i,'chrom']

        cmd = '{0} {1} {2} {3} {4} {5}'.format('bigWigSummary',bw, chrom, st, end, end-st)

        result = os.popen(cmd)
        content = result.readlines()
        temp = content[0].strip().replace('n/a','0')
        temp = temp.split('\t')
        values = np.array(temp,dtype=np.float64)

        values  = np.hstack((np.zeros(st - cnt + pad),values,np.zeros(cnt + pad - end +1 )))

        df.loc[i,'RP'] = np.dot( values, wt)
    
    return df


dp = bigWig2bed_pot(bw, bd,'mm9',pad=1000)

print(dp.head(30))
