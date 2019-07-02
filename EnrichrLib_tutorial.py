
# coding: utf-8

# In[1]:

import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

# # Import ErichRLib

# In[2]:

import EnrichRLib as erl

def sigmoid(x):                                        
   return 1 / (1 + np.exp(-x))


# Gene list may be imported from file or direct entered in code as list object

# In[3]:

gl = [  'CD3E', 'BLK', 'PTPN22', 'PAG1', 'CTLA4', 'PIK3CD', 'LAT2', 'CSK', 
        'CD247', 'CD3G', 'THEMIS', 'PSMB8', 'LCP2', 'GATA3', 'LAT', 'SLA2', 
        'SKAP1', 'TRAT1', 'BCL2', 'CD3D', 'THY1','RUNX1', 'BLK', 'PTPN22']


# # Fetching gene sets from Enrichr website 

# In[4]:

# sub directory where gene sets in .gmt format will be placed
lib_dir='_tmp'

# List of gene sets to fetch, by default all gene sets will be fetched
gss = [ 
       'GO_Biological_Process_2018',
       'GO_Cellular_Component_2018',
       'GO_Molecular_Function_2018',
       'KEGG_2016',
       'Reactome_2016'
       ]

#~ erl.get_Enrichr(out_dir='_tmp', libs=gss)


# # Enrichment of single gene set

# In[5]:

# gene set name and file name
#~ lib_dir='_tmp'
#~ gsn = 'Reactome_2016'
#~ gs_fn =lib_dir+'/'+gsn+'.gmt'

#~ # read gmt to python dictionary
#~ gs = erl.read_gmt(gs_fn)

#~ print(gs['Formation of Incision Complex in GG-NER_Homo sapiens_R-HSA-5696395'])


# ### Enrichment run

# In[6]:

#~ enr = erl.enrich(gl, gs)


# In[7]:

# Plot result as barplot and compare with EnrichR image

import seaborn as sns

#~ ds = enr.head(10)
#~ sns.barplot(y=ds.index,
                #~ x='-log10(p-Val)', 
                #~ color="Red",
                #~ data = ds)
#~ plt.title(gsn)


# ![title](im/React2016_EnrichTest_bar.png)

# In[8]:

#~ ds


# ![title](im/React2016_EnrichTest_table.png)

# # Enrichment of several gene sets: batch analysis

# In[9]:

# List of gene sets as above
#~ gss = [ 
       #~ 'GO_Biological_Process_2018',
       #~ 'GO_Cellular_Component_2018',
       #~ 'GO_Molecular_Function_2018',
       #~ 'KEGG_2016',
       #~ 'Reactome_2016'
       #~ ]


# ### Batch enrichment 

# In[10]:

enrr = erl.enrich_gs(gl,gss, path_lib=lib_dir)


# ### Plots 

# In[11]:

#~ enrr.sort_values('p-Val', axis=0, inplace = True)
#~ ds = enrr.head(20)

#~ f, ax = plt.subplots()
#~ sns.barplot(y=ds.index,
            #~ x='-log10(p-Val)',
            #~ ax = ax, 
            #~ color="Red", 
            #~ data = ds)
#~ ax.set_title('All terms')


# # Clustering 

# In[12]:

# For futher analysis it is convinient to filter terms by p-value for
enrr = enrr[enrr['p-Val']<0.0005]
#~ len(enrr)


# ## Calculate closeness by kappa-score

# In[13]:

#~ dd = erl.kappa_matrx(gl,enrr)


# In[14]:

## Show cluster based on kappa-score closeness

#~ grid = sns.clustermap(dd, cmap='Blues', figsize=(16, 16))


# ## Clastering:  top level

# In[15]:

## Cluster: this calculate and add cluster number column
#~ enrr = erl.cluster(gl,enrr)

# Make additional calculation on existing columns for visualization
#~ enrr.loc[:,'ass_genes_percnt'] = 100*enrr.loc[:,'num_list']/enrr.loc[:,'num_term']
#~ enrr.sort_values('cluster', axis=0, inplace = True)

# use consistent discrete palette
cm = 'tab20'

#~ ds = enrr.head(40)

#~ f, ax = plt.subplots(figsize=(18.5, 8.5))
#~ sns.barplot(y=ds.index,
            #~ x='-log10(p-Val)',
            #~ ax = ax, 
            #~ hue ='cluster',
            #~ dodge=False,
            #~ data = ds,
            #~ palette = cm)
#~ ax.set_title('All terms')


# In[16]:

#~ f, ax = plt.subplots(figsize=(18.5, 8.5))
#~ sns.barplot(y=ds.index,
            #~ x='ass_genes_percnt',
            #~ ax = ax, 
            #~ hue ='cluster', 
            #~ dodge=False,
            #~ data = ds,
            #~ palette = cm)
#~ ax.set_title('All terms')
#~ ax.set_xlabel('%Genes/Term')


# # Network construction based on cluster

# In[17]:

## Graph and network
## control conectivity by kappa score parameter, default kappa=0.4
G, enr_out, nt = erl.make_graph(gl, enrr, kappa=0.8)

#~ spring = 100
#~ for u,v in G.edges:
    #~ G.edges[u, v]['spring'] = sigmoid(G.edges[u, v]['weight'])*(spring/500)

#~ pos=nx.spring_layout(G,k=0.1, weight = 'spring')

#~ df = pd.DataFrame.from_dict(pos, orient='index',columns=['x','y'])
#~ cl = [G.nodes[v]['cluster'] for v in df.index]
#~ sz = [G.nodes[v]['mlog10pVal'] for v in df.index]


#~ sns.scatterplot(x='x',
                #~ y='y',
                #~ data=df, 
                #~ hue=cl,
                #~ legend  ="full",
                #~ size = sz, 
                #~ sizes=(20,400),
                #~ hue_norm =(1,20),
                #~ palette=cm)


erl.draw_graph(G, spring=5, pval_prcnt=0.7, palette= plt.cm.tab20)



enr_out.sort_values('cluster', axis=0, inplace = True)

# use consistent discrete palette


dss = enr_out.copy()#head(40)

f, ax = plt.subplots(figsize=(18.5, 8.5))
sns.barplot(y=dss.index,
            x='-log10(p-Val)',
            ax = ax, 
            hue ='cluster',
            dodge=False,
            data = dss,
            palette = cm)
ax.set_title('All terms')



plt.show()
# In[18]:

#~ #enrr.head(40).index
#~ enrr.loc['T cell receptor signaling pathway (GO:0050852)']


#~ # In[29]:

#~ G.nodes['T cell receptor signaling pathway (GO:0050852)']['cluster']


# In[ ]:




# ## Exporting network tables for Cytoscape

# In[19]:

# If save result and network for Cytoscape
#~ exp_name='GO_KEGG_React'

#~ SAVE = False

#~ if SAVE:
    #~ enrr.to_csv(exp_name+'-enrich.txt', sep="\t")
    #~ nt.to_csv(exp_name+'-network.txt', sep="\t")
    #~ nt_tb.to_csv(exp_name+'-network_table.txt', sep="\t")


# ![title](im/React2016_EnrichTest_net.png)

# # ClueGO for same gene list and coresponding gene sets

# ![title](im/React2016_EnrichTest_netClue.png)

# In[68]:




# In[69]:

#lbm
#for v in G:
#    print(G.nodes[v]['num_list'])


# In[79]:

#~ import networkx as nx
#~ cl = [G.nodes[v]['cluster'] for v in G] 
#~ sz = [G.nodes[v]['mlog10pVal']*100 for v in G] 
#~ shr=10
# lb_M = {v:v for v in G if G.nodes[v]['num_list']/G.nodes[v]['num_term']>shr}
# lb_m = {v:v for v in G if G.nodes[v]['num_list']/G.nodes[v]['num_term']<=shr}
#~ lb_M = {v:v for v in G if G.nodes[v]['mlog10pVal']>shr}
#~ lb_m = {v:v for v in G if G.nodes[v]['mlog10pVal']<=shr}



#~ pos=nx.spring_layout(G,k=.2, iterations=100)
#~ f1, ax1 = plt.subplots(figsize=(20, 20))
   #~ #pos=nx.kamada_kawai_layout(G)  
#~ nx.draw(G, pos=pos,
        #~ #with_labels = True,
        #~ labels = lb_M,
        #~ node_color = cl,
        #~ node_size =  sz,
        #~ font_size=14,
        #~ font_weight='bold',
        #~ cmap = cm) 


#~ nx.draw(G, pos=pos,
        #~ #with_labels = True,
        #~ labels = lb_m,
        #~ node_color = cl,
        #~ node_size =  sz,
        #~ font_size=6,
        #~ cmap = cm)


# In[38]:

#~ import numpy as np
#~ np.argmax(sz)


# In[39]:

#~ print(cl[57])


# In[ ]:



