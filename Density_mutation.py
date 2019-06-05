import numpy as np
import scipy
import pandas as pd
import pybedtools as pb


import matplotlib.pyplot as plt
import seaborn as sns

from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')
WGS = join(DATADIR,'tracks/WGS-WES/Germline')


import allel
import deeptools.getScorePerBigWigBin as gs
from scipy.signal import savgol_filter


# my libs
from Enh_Mut_Manip import variants_bed_counts

#=================

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




def bed_bins(step, genome='mm9'):
    gnm = pb.chromsizes(genome)
    all_bins = pd.DataFrame()

    for chrm  in gnm.keys():
        st = np.arange(gnm[chrm][0],gnm[chrm][1],step)
        bins = pd.DataFrame()
        bins['start'] = st[:-1]
        bins['end'] = st[1:]
        bins['chrom'] = chrm
        bins['score'] = 0
        all_bins = pd.concat([all_bins, bins])
        
    cols = ['chrom','start','end','score']

    all_bins = all_bins[cols]
    bin_bed = pb.BedTool.from_dataframe(all_bins)

    return bin_bed




def hypmut_bw(vcf,bw,step=1000,nstep=50,shr=0.8,genome='mm9'):
    gnm = pb.chromsizes(genome)
    bin_bed = bed_bins(step, genome=genome)

    vcf_b = pb.BedTool(vcf)
    bed_out = variants_bed_counts(vcf_b,bin_bed)
    out = bed_out.to_dataframe()
    out['density'] = out['score']/step

    hyp_bw = pd.DataFrame()
    
    for chrom  in gnm.keys():
        out_chr = out[out['chrom']==chrom]
        x = out_chr['start'].values  # + int(step/2)
        y = out_chr['density'].values
        
        st = gnm[chrom][0]
        end = gnm[chrom][1]
        
        vl = gs.countFragmentsInRegions_worker(chrom, 
                                                int(st), 
                                                int(end), 
                                                [bw], 
                                                stepSize=step,
                                                binLength=step, 
                                                save_data=False)


        xb = np.arange(st,end,step)[:-1]
        yb = np.squeeze(vl[0])[:-1]
        
        x_r = x[:-nstep]

        cr_f = np.zeros(x_r.shape)

        for i in range(len(x_r)):
            cr_f[i] = np.corrcoef(y[i:i+nstep],yb[i:i+nstep])[0,1]

        x_c = x_r[cr_f>shr]

        hyp_cor= pd.DataFrame()
        hyp_cor['start']= x_c
        hyp_cor['end']= x_c+nstep*step
        hyp_cor['chrom']= chrom 

        hyp_cor = hyp_cor[['chrom','start','end']].astype({'start':int,'end':int})
        hyp_bw = pd.concat([hyp_bw, hyp_cor])
    
    hyp_bed = pb.BedTool.from_dataframe(hyp_bw).merge()
    return hyp_bed




#~ chrom = 'chr1'
#~ step = 2000

# VCF
vcf = pb.BedTool(join(WGS,'TLX3_WGS.vcf.gz'))

# BigWig track
trck_list = [
        'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_H3K4me1_mm9.bw',
        'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_H3K4me2_mm9.bw',
        'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_H3K4me3_mm9.bw',
        'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_H3K9ac_mm9.bw',
        'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_H3K9me3_mm9.bw',
        'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_H3K27ac_mm9.bw',
        'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_H3K27me3_mm9.bw',
        'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_H3K36me3_mm9.bw',
        'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_POLII_mm9.bw',
        'tracks/ChiP-seq_tracks/mm9_bigWig/TLX3_TLX3_mm9.bw'
]


#~ for trck in trck_list: 
    #~ bw = join(DATADIR,trck)                

    #~ hyp_bdd = hypmut_bw(vcf,bw, step=600, nstep=150,shr=0.85)
    #~ hyp_bdd.saveas(join(WGS,'Hypermut_TLX3_WGS_'+trck.split('/')[-1]+'.bed'))




res_trck_l = [
        'Hypermut_TLX3_WGS_TLX3_H3K4me1_mm9.bw.bed',
        'Hypermut_TLX3_WGS_TLX3_H3K4me2_mm9.bw.bed',
        'Hypermut_TLX3_WGS_TLX3_H3K4me3_mm9.bw.bed',
        'Hypermut_TLX3_WGS_TLX3_H3K9ac_mm9.bw.bed',
        'Hypermut_TLX3_WGS_TLX3_H3K9me3_mm9.bw.bed',
        'Hypermut_TLX3_WGS_TLX3_H3K27ac_mm9.bw.bed',
        'Hypermut_TLX3_WGS_TLX3_H3K27me3_mm9.bw.bed',
        'Hypermut_TLX3_WGS_TLX3_H3K36me3_mm9.bw.bed',
        'Hypermut_TLX3_WGS_TLX3_POLII_mm9.bw.bed',
        'Hypermut_TLX3_WGS_TLX3_TLX3_mm9.bw.bed'
]


allb = pb.BedTool(join(WGS,res_trck_l[0]))

for trck in res_trck_l[1:]: 
    btr = pb.BedTool(join(WGS,trck))
    allb = allb+btr
    



#~ hyp_cor= pd.DataFrame()

# Figure
step=2000

bin_bed = bed_bins(step)





bed_out = variants_bed_counts(vcf,bin_bed)


out = bed_out.to_dataframe()
out['density'] = out['score']/step


chrom='chr1'

out_chr = out[out['chrom']==chrom]

x = out_chr['start'].values + step/2
y = out_chr['density'].values

#yv = savgol_filter(y,13,5)
yv = y.copy()


trck = trck_list[5]

bw = join(DATADIR,trck)

mm9 = pb.chromsizes('mm9')

st = mm9[chrom][0]
end = mm9[chrom][1]

vl = gs.countFragmentsInRegions_worker(chrom, 
                                        int(st), 
                                        int(end), 
                                        [bw], 
                                        stepSize=step,
                                        binLength=step, 
                                        save_data=False)


xb = np.arange(st,end,step)[:-1]
yb = np.squeeze(vl[0])[:-1]


#~ # Figures
f, ax = plt.subplots(figsize=(26.5, 6.5))

ax.plot(x,yv,'--.', color='b', MarkerSize=0.8, LineWidth=0.1)
ax.set_ylabel('Hypermutation density')

ax1 = ax.twinx()
ax1.plot(xb,yb, '--.',color='r', MarkerSize=0.8, LineWidth=0.1)

ax1.set_ylabel(trck.split('/')[-1])


#~ corr  = np.corrcoef(yv,yb)

#~ print(corr)

#~ f2, ax2 = plt.subplots()
#~ ax2.scatter(yv,yb, s=1)




nstep = 50
x_r = x[:-nstep]

cr_f = np.zeros(x_r.shape)

for i in range(len(x_r)):
    cr_f[i] = np.corrcoef(yv[i:i+nstep],yb[i:i+nstep])[0,1]


f3, ax3 = plt.subplots()
ax3.plot(x_r,cr_f, '--.',color='r', MarkerSize=0.8, LineWidth=0.1)

plt.show()

# Pocessing

#~ x_c = x_r[cr_f>0.8]


#~ hyp_cor= pd.DataFrame()

#~ hyp_cor['start']= x_c

#~ hyp_cor['end']= x_c+nstep*step

#~ hyp_cor['chrom']= chrom 

#~ hyp_cor = hyp_cor[['chrom','start','end']].astype({'start':int,'end':int})

#~ hyp_bed = pb.BedTool.from_dataframe(hyp_cor).merge()


#~ plt.show()
