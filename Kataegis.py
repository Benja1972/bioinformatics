#Kataegis


import numpy as np
import scipy
import pandas as pd
import pybedtools as pb


#import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')
WGS = join(DATADIR,'tracks/WGS-WES/Germline')


# Load table
import allel


#vcf_tb = allel.vcf_to_dataframe(join(WGS,'Ehn_Active_TLX3_mut.vcf'),fields='*', numbers={'ALT': 1}, transformers=allel.ANNTransformer())
vcf_tb = allel.vcf_to_dataframe(join(WGS,'TLX3_WGS.vcf.gz'),fields='*', numbers={'ALT': 1}, transformers=allel.ANNTransformer())

cols = ['CHROM', 'POS', 'REF', 'ALT','is_snp']

vcf_tb = vcf_tb[cols]

def kataegis(chrm,vcf_tb):

    vcf_ch = vcf_tb[vcf_tb['CHROM']==chrm].sort_values('POS',axis=0, ascending=True)

    vcf_ch['DIFF'] = vcf_ch['POS'].diff()
    vcf_ch['logDIFF'] = np.log2(vcf_ch['DIFF'])
    vcf_ch['mut_TYPE'] = vcf_ch['REF']+'>'+vcf_ch['ALT']
    vcf_ch['mut_TYPE'][vcf_ch['is_snp']==False]='not_snp'

    f, ax = plt.subplots(figsize=(16.5, 6.5))
    hue_order  = ['A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G', 'not_snp']
    sns.scatterplot(x='POS', 
                    y = 'DIFF', 
                    hue='mut_TYPE',
                    hue_order = hue_order,
                    palette= 'tab20_r', #gnuplot_r', #tab20',
                    ax = ax,
                    data=vcf_ch)

    ax.set_title(chrm)


for chrm in  list(vcf_tb['CHROM'].unique())[:2]:
    kataegis(chrm,vcf_tb)



#f, ax = plt.subplots(2, 1, sharex=True)
#f.subplots_adjust(hspace=0)

#~ ax.eventplot(vcf_1['POS'], 
                #~ colors='k', 
                #~ linelengths=1000,
                #~ lineoffsets=100)


#~ ax.scatter(x = vcf_1['POS'], 
            #~ y = vcf_1['DIFF']) 
            #~ #c=vcf_1['mutTYPE'], 
            #~ #label=vcf_1['mutTYPE'])


#ax[0].spines['bottom'].set_visible(False)
#ax[0].tick_params(labeltop=False)



#ax[1].set_axis_off()
#~ ax[1].axis['top'].set_visible(False)
#~ ax[1].axis['left'].set_visible(False)
#~ ax[1].axis['right'].set_visible(False)
#ax[1].set_autoscale_on(False)
#~ #ax[1].set_axisbelow(False)
#~ ax[1].set_ylim(0, 0.1)


plt.show()
