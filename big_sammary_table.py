
import numpy as np
import yaml
import os
from os.path import join 
import pybedtools as pb
import gffutils
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import matplotlib.pyplot as plt
import metaseq
import pandas as pd


ref_dir = '../../references/mm9'
data_dir = 'tracks'
rel_path = "/home/sergio/media"



file_t = open(join(data_dir,"TLX3_list.yaml"), "r")
file_r = open(join(data_dir,"RAG_list.yaml"), "r")

tlx_lst = yaml.load(file_t)
rag_lst = yaml.load(file_r)


tsses     = pb.BedTool(join(ref_dir,'mm9_tsses.gtf'))
tsses_1kb = pb.BedTool(join(ref_dir,'mm9_tsses-1kb.gtf'))
tsses_3kb = pb.BedTool(join(ref_dir,'mm9_tsses-3kb.gtf'))
tsses_5kb = pb.BedTool(join(ref_dir,'mm9_tsses-5kb.gtf'))

# ChiP-seq signals around TSS loading
features, arrays = metaseq.persistence.load_features_and_arrays(prefix='tlx-list1K')
features3K, arrays3K = metaseq.persistence.load_features_and_arrays(prefix='tlx-list3K')
features_rag, arrays_rag = metaseq.persistence.load_features_and_arrays(prefix='rag-list1K')
features3K_rag, arrays3K_rag = metaseq.persistence.load_features_and_arrays(prefix='rag-list3K')



from metaseq.results_table import ResultsTable, DESeq2Results

# RNA-seq expression tables
tlx_tabf = rel_path+tlx_lst["table"][0]
tlx_vs_rag_table = ResultsTable(tlx_tabf, import_kwargs=dict(index_col=0))

# Add gene names for first column
gene_names = pd.read_table("tracks/UCSC_mm9_transcripID_to_geneSymbol.sort.txt", index_col=0,names=["Geneid", "Gene_name"])

tlx_vs_rag_table.data = gene_names.join(tlx_vs_rag_table.data, how='right')

# Mean values as new columns
tlx_vs_rag_table.data['mean(TLX3)'] = (tlx_vs_rag_table.data['TLX3.1_1']+tlx_vs_rag_table.data['TLX3.1_5']+tlx_vs_rag_table.data['TLX3.1_P'])/3.
tlx_vs_rag_table.data['mean(RAG)'] = (tlx_vs_rag_table.data['R2.RAG1W.RAG1']+tlx_vs_rag_table.data['RAGS.RAGZ']+tlx_vs_rag_table.data['RAGZ'])/3.


# Re-align the ResultsTables to match the GTF file
tlx_vs_rag_table = tlx_vs_rag_table.reindex_to(tsses, attribute='transcript_id')

old = [
    'rag-POLII',
    'rag-H3K27me3',
    'rag-H3K4me2',
    'rag-H3K4me3',
    'rag-H3K9ac',
    'rag-H3K36me3',
    'rag-H3K9me3']


# ChiP-Seq signal in charge
for chname in  arrays.keys():
    chnameR = 'rag'+chname[3:]

    charr = arrays[chname]
    charrR = arrays_rag[chnameR]
    charr3 = arrays3K[chname]
    charrR3 = arrays3K_rag[chnameR]
    
    if any(chnameR in string for string in old):
        k = (abs(np.max(charr3)-np.min(charr3)))/(abs(np.max(charrR3)-np.min(charrR3)))
        charrR = charrR*k
        charrR3 = charrR3*k
        chnameR = chnameR+'_old'
    
    tlx_vs_rag_table.data[chname+'_TSS_mean1K'] = charr.mean(axis=1)
    tlx_vs_rag_table.data[chname+'_TSS_mean3K'] = charr3.mean(axis=1)
    
    tlx_vs_rag_table.data[chnameR+'_TSS_mean1K'] = charrR.mean(axis=1)
    tlx_vs_rag_table.data[chnameR+'_TSS_mean3K'] = charrR3.mean(axis=1)


## Save to file
# tlx_vs_rag_table.data.to_csv('test.csv')







# # Top 200 from ChiP-seq


## Sort features for feature analysis. Not used for below
# sorting features by TLX strength around TSS
#~ tlx_mean = charr.mean(axis=1)
#~ tlx_mean_srt = np.sort(tlx_mean)[::-1]
#~ sort_indices = np.argsort(tlx_mean)[::-1]

#gene_names = pd.read_table("tracks/UCSC_mm9_transcripID_to_geneSymbol.sort.txt", index_col=0,names=["Geneid", "Gene_name"])

#gene_names.head(20)


#top200names = gene_names.loc[top200.index]



#top200all = top200names.join(top200)











#======================================================================
#======================================================================
#======================================================================
    


