#!/usr/bin/env python3

import pybedtools as pb
import numpy as np
import yaml
from os.path import join 
import pandas as pd

import matplotlib.pyplot as plt

SAVE = False

## === Bed manipulation
enh= pb.BedTool('tracks/TLX3_6_FE_E4_sorted.bed')
enh_df = enh.to_dataframe()
enh_df['len'] = enh_df['end'] - enh_df['start']
enh_df['enh_name'] = enh_df.index
enh_df['enh_name'] = 'enh_'+enh_df['enh_name'].astype(str)
enh_df.drop('name', axis=1, inplace=True)


enh_sh = enh_df[enh_df['len']<=4500]
enh_ln = enh_df[enh_df['len']>4500]


plt.hist(enh_df['len'], bins=50)


col2bed = ['chrom', 'start', 'end', 'enh_name']


if SAVE:
    enh_bd =  pb.BedTool.from_dataframe(enh_df[col2bed])
    enh_bs =  pb.BedTool.from_dataframe(enh_sh[col2bed])
    enh_bl =  pb.BedTool.from_dataframe(enh_ln[col2bed])
    
    enh_bd.saveas('tracks/enhancers_ChromHMM_all.bed')
    enh_bs.saveas('tracks/enhancers_ChromHMM_short.bed')
    enh_bl.saveas('tracks/enhancers_ChromHMM_long.bed')


plt.show()



# Enhancers Table manipulation
gn2enh = pd.read_table('tracks/enhancers_ChromHMM-gene2enh.txt', 
                        header=1, 
                        names=['gene_name','enhancers'])


gl=['Dusp4', 'Ms4a6b' ,'Rgs18','Tha1','Aip']


df = pd.DataFrame()


for gen in gl:
    reg = gn2enh[gn2enh['gene_name']==gen]

    dr = reg.iloc[0]['enhancers'].split(' ')
    enhl = [x for x in dr if 'enh' in x]

    ehGene = enh_df.loc[enh_df['enh_name'].isin(enhl)]
    ehGene['enh_name'] = ehGene["enh_name"].map(str) + '|' + gen
    df = pd.concat([df,ehGene], ignore_index=True)


#enhm_ = pb.BedTool('tracks/enhancers_ChromHMM_all.bed')
#enhDf = enhm_Nkb.to_dataframe()


