import numpy as np
from os.path import join 
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import EnrichRLib as erl
from matplotlib.backends.backend_pdf import PdfPages


# === User defined
import RNA_expression_processing as rn


def log2p1(x):
    return np.log2(x + 1)



SAVE = True


# === Load file

path_out = 'tracks/MARGE/relativeRP/result/'



### STOP ###
#~ plt.show()
#~ sys.exit(0)


A = 'TLX3'
B = 'RAG'

## Up analysis
fup = path_out+'UP_TLX3_peaks'+'_'+A+'_vs_'+B+'.csv'
fep = path_out+'UP_RP_enh_TLX3_pck_genes.txt'

df_fup = pd.read_csv(fup, index_col=0)
gup = list(df_fup['name'].str.upper())

df_fep = pd.read_table(fep, names=['name'])
gep = list(df_fep['name'].str.upper())


guep = list(set(gup)|set(gep))





#~ if SAVE:
    #~ dn_rp_gene.to_csv(path_out+'DN_TLX3_peaks'+'_'+A+'_vs_'+B+'.csv')

### =================





### =================


### == DN analysis

    #~ with open(path_out+'DN_RP_enh_TLX3_pck_genes.txt', 'w') as fp:
        #~ fp.write("\n".join(dn_enh_tlx_genes))


### =================











### === Enrichr ======
### ==================

gss = [
    'BioCarta_2016',
    'Cancer_Cell_Line_Encyclopedia',
    'ChEA_2016',
    'Disease_Perturbations_from_GEO_down',
    'Disease_Perturbations_from_GEO_up',
    'Disease_Signatures_from_GEO_down_2014',
    'Disease_Signatures_from_GEO_up_2014',
    'GO_Biological_Process_2018',
    'GO_Cellular_Component_2018',
    'GO_Molecular_Function_2018',
    'Human_Phenotype_Ontology',
    'KEGG_2016',
    'MSigDB_Computational',
    'MSigDB_Oncogenic_Signatures',
    'Mouse_Gene_Atlas',
    'NCI-60_Cancer_Cell_Lines',
    'NCI-Nature_2016',
    #'OMIM_Disease',
    'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO',
    'Reactome_2016',
    'WikiPathways_2016'
        ]

    




enrDD  = erl.enrich_gs(gup,gss)

enrDD.sort_values('p-Val', axis=0, inplace = True)
ds = enrDD.head(20)
f, ax = plt.subplots(figsize=(16.5, 6.5))
sns.barplot(y=ds.index,
            x='-log10(p-Val)',
            ax = ax, 
            color="Red", 
            data = ds)
ax.set_title('All terms')

plt.show()







