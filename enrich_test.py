from os.path import join
import numpy as np
import pandas as pd
#import RNA_expression_processing as rn
import EnrichRLib as erl
import seaborn as sns
import matplotlib.pyplot as plt


SAVE = False

exp_name='GO_KEGG2'


gt = pd.read_table('tracks/MARGE/relativeRP/DN_RP_enhancers_genes.txt')


setR = set(list(gt.iloc[:,0]))







gss = [ 
        #~ 'Pierre_gene_sets.gmt',
        'GO_Biological_Process_2018.gmt',
       'GO_Cellular_Component_2018.gmt',
       'GO_Molecular_Function_2018.gmt',
       'KEGG_2016.gmt',
       'Reactome_2016.gmt'
       ]



enrr = pd.DataFrame()



for gs in gss:
    pl = erl.read_gmt('EnrichrLibs/'+gs)
    enr = erl.enrich(setR, pl)
    print(len(enr))
    enrr = pd.concat([enrr, enr])
    #~ ds = enr.head(20)
    #~ f, ax = plt.subplots(figsize=(16.5, 6.5))
    #~ sns.barplot(y=ds.index,
                #~ x='-log10(p-Val)',
                #~ ax = ax, 
                #~ color="Red", 
                #~ data = ds)
    #~ ax.set_title(gs.split(".")[0])




# p-Value filter
enrr = enrr[enrr['p-Val']<0.005]

## Show cluster
dd = erl.dist_matrx(setR,enrr)
grid = sns.clustermap(dd, cmap='Blues', figsize=(16, 16))



enrr = erl.cluster(setR,enrr)


if SAVE:
    enrr.to_csv(exp_name+'-enrich.txt', sep="\t")


nt, nt_tb, G = erl.make_graph(setR, enrr, draw=True)


if SAVE:
    nt.to_csv(exp_name+'-network.txt', sep="\t")
    nt_tb.to_csv(exp_name+'-network_table.txt', sep="\t")


plt.show()
