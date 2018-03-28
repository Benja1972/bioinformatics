## Minimal test of GSEAPY
## Dev version of GSEAPY, run in Python3 
## $ conda install -c bioconda gseapy
## conda env "gseapy_test --> /home/sergio/miniconda3/envs/gseapy_test"


import gseapy as gp
import numpy as np
from os.path import join 
import pandas as pd
import matplotlib.pyplot as plt



tbl = pd.read_table(join('tracks', 'TLX3vsRAG-results_genes.txt'), index_col=0)

tbl = tbl[(tbl.padj < 0.05)].dropna()
names = pd.read_table("tracks/annot_tracks/references/mm9/mm9_EnsemblTransc_GeneNames.txt", 
                           index_col=0,
                           header=0,
                           names=['GeneID', 'TransID', 'Gene_name'])

names = names.drop('TransID', axis=1).drop_duplicates()
names = names.reindex(tbl.index)


cols = ['TLX3.1_1', 'TLX3.1_5', 'TLX3.1_P', 'R2.RAG1W.RAG1', 'RAGS.RAGZ', 'RAGZ']
classes = ['TLX3','TLX3','TLX3', 'RAG','RAG','RAG']

tbn = pd.concat([names,tbl[cols]], axis=1)

#tbn = tbn.drop(['baseMean', 'log2FoldChange', 'lfcSE', 'stat',  'pvalue',  'padj'], axis=1)

tbn['Gene_name'] = tbn['Gene_name'].str.upper()
tbn = tbn.reset_index(drop=True)



## ==== GSEAPY
# --- get gene lists 
#~ list_gs = gp.get_library_name()

gs_dic = {  'Hallmark': 'tracks/GSEA_gene_sets/h.all.v6.0.symbols.gmt',
            'BioCarta_2013': 'BioCarta_2013',
            'BioCarta_2015': 'BioCarta_2015',
            'BioCarta_2016': 'BioCarta_2016',
            'KEGG_2013': 'KEGG_2013',
            'KEGG_2015': 'KEGG_2015',
            'KEGG_2016': 'KEGG_2016',
            'Reactome_2013': 'Reactome_2013',
            'Reactome_2015': 'Reactome_2015',
            'Reactome_2016': 'Reactome_2016',
            'Canonical_pathways': 'tracks/GSEA_gene_sets/c2.cp.v6.0.symbols.gmt',
            'Imunno':'tracks/GSEA_gene_sets/c7.all.v6.0.symbols.gmt',
            'Oncogenic_signatures': 'MSigDB_Oncogenic_Signatures',
            'Computational':'MSigDB_Computational',
            'NCI-60_Cancer_Cell_Lines': 'NCI-60_Cancer_Cell_Lines',
            'Transcription_factor_targets': 'tracks/GSEA_gene_sets/c3.tft.v6.0.symbols.gmt'}




#~ for g_set in gs_dic.keys():
for g_set in ['Transcription_factor_targets']:
    out_dir = 'GSEA/TLX3vsRAG_' + g_set

    gs_res = gp.gsea(data=tbn, 
                    gene_sets = gs_dic[g_set], 
                    outdir=out_dir, 
                    cls = classes)
    
    # plotting
    gsea_results= gs_res.res2d
    with plt.style.context('ggplot'):
        gsea_results = gsea_results.reset_index()
        gsea_results.head(40).plot.barh(y='fdr',x='Term', figsize=(18, 6), fontsize=10)
        plt.gca().invert_yaxis()



    plt.savefig(out_dir+'/'+'TLX3vsRAG_' +g_set+'.pdf', format="pdf")


plt.show()
