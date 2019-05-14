import numpy as np
import scipy
import pandas as pd
import pybedtools as pb


#import matplotlib as mpl
#~ import matplotlib
#~ matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns

from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')
WGS = join(DATADIR,'tracks/WGS-WES/Germline')


def kataegis(chrm,vcf_tb, shr=0):

    vcf_ch = vcf_tb[vcf_tb['CHROM']==chrm].sort_values('POS',axis=0, ascending=True)
    
    vcf_ch['mut_TYPE'] = vcf_ch['REF']+'>'+vcf_ch['ALT']
    vcf_ch['mut_TYPE'][vcf_ch['is_snp']==False]='not_snp'
    
    vcf_ch['DIFF'] = vcf_ch['POS'].diff()
    vcf_ch['logDIFF'] = np.log10(vcf_ch['DIFF'])
    if shr > 0:
        vcf_ch = vcf_ch[vcf_ch['logDIFF']<shr]

    f, ax = plt.subplots(figsize=(26.5, 6.5))
    hue_order  = ['A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G', 'not_snp']
    sns.scatterplot(x= 'POS', #range(len(vcf_ch)), 
                    y = 'logDIFF', 
                    hue='mut_TYPE',
                    hue_order = hue_order,
                    palette= 'tab20_r', #gnuplot_r', #tab20',
                    ax = ax,
                    data=vcf_ch)

    ax.set_title(chrm)
    
    return f, ax



chrm = 'chr1'

# WGS track
import allel

vcf_tb = allel.vcf_to_dataframe(join(WGS,'TLX3_WGS.vcf.gz'),fields='*', numbers={'ALT': 1})

vcf_tb = vcf_tb[vcf_tb['FILTER_PASS']==True]

cols = ['CHROM', 'POS', 'REF', 'ALT','is_snp']

vcf_tb = vcf_tb[cols]




#vcf_tb = vcf_tb[(vcf_tb['CHROM']==chrm) & (vcf_tb['POS']<x.max()) & (vcf_tb['POS']>x.min())]
#vcf_tb = vcf_tb[(vcf_tb['CHROM']==chrm) & (vcf_tb['POS']<x.max()) & (vcf_tb['POS']>x.min())]

f, ax = kataegis(chrm,vcf_tb)


# test
import deeptools.getScorePerBigWigBin as gs

bw = join(DATADIR,'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_H3K27ac_mm9.bw')



#~ chrm = 'chr1'
st = 0
end = vcf_tb['POS'].max() #1.5e8
step=200

vl = gs.countFragmentsInRegions_worker(chrm, 
                                        int(st), 
                                        int(end), 
                                        [bw], 
                                        stepSize=step,
                                        binLength=step, 
                                        save_data=False)


xv = np.arange(st,end,step)
yv = np.squeeze(vl[0])
yvf = savgol_filter(yv,10,3)



ax1 = ax.twinx()
colr = 'r'
ax1.plot(xv,yvf, '--.',color=colr, MarkerSize=0.8, LineWidth=0.1)
ax1.set_ylim(0, 3)
ax1.tick_params(axis='y', labelcolor=colr)



plt.show()
