import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl



bd  = '/home/sergio/Res_CIML/TLX3_project/data/tracks/ChromHMM_enhancers/ChromHMM03_state_6/RAG_6_dense.bed'

st = pd.read_table(bd, header=1, names=['chr', 'start', 'end', 'state', 'score', 'strand','str', 'ed','color' ])


colr = st[['state', 'color']]
colr.drop_duplicates(inplace=True)
colr.sort_values('state', inplace=True)

cmap = mpl.colors.ListedColormap(np.array(list(colr['color'].str.split(','))).astype(int)/256)




fig = plt.figure(figsize=(3,len(colr)))

ax1 = fig.add_axes([ 0.1, 0.15, 0.15, 0.80])


bounds = list(colr['state'])+[len(colr)+1]
mpl.colorbar.ColorbarBase(ax1, 
                        cmap=cmap, 
                        boundaries= bounds,
                        orientation='vertical')

plt.show()


