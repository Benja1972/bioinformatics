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


def closest_genes(bed,genes):
    bed = bed.sort()
    near = bed.closest(genes, D='ref')
    near_d = near.to_dataframe()
    
    gl = list(near_d['thickStart'].str.upper().unique())
    
    return gl




bd = pb.BedTool(join(DATADIR,'tracks/WGS-WES/Germline/Hypermut_TLX3_WGS_TLX3_TLX3_mm9.bw.bed'))
tss = pb.BedTool(join(DATADIR,'tracks/annot_tracks/references/mm9/mm9.refGene.bed'))    # [1]



gl =  closest_genes(bd,tss)

import EnrichRLib as erl

## Enrichment analysis
# List of gene sets 
gss = [ 
        #~ 'GO_Biological_Process_2018',
        #~ 'GO_Cellular_Component_2018',
        #~ 'GO_Molecular_Function_2018',
        #~ 'KEGG_2016',
        #~ 'Reactome_2016',
        'Cancer_Cell_Line_Encyclopedia',
        'NCI-60_Cancer_Cell_Lines',
        'MSigDB_Computational',
        'MSigDB_Oncogenic_Signatures'
       ]





enr = erl.enrich_gs(gl,gss, path_lib='../data/EnrichrLibs')




# For futher analysis it is convinient to filter terms by p-value
enr_a = enr[enr['p-Val']<0.05]

G, enr_out, nt = erl.make_graph(gl, enr_a, kappa=0.4)


enr_out.sort_values('cluster', axis=0, inplace = True)

# --- Plot ---
cm = 'tab20'
erl.draw_graph(G, spring=500, pval_prcnt=0.7, palette= cm)


ds = enr_out.head(40)

f, ax = plt.subplots(figsize=(12, 12))
sns.barplot(y=ds.index,
            x='-log10(p-Val)',
            ax = ax, 
            hue ='cluster',
            dodge=False,
            data = ds,
            palette = cm)
ax.set_title('TLX3 peaks with hypermutations')

plt.show()
