#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
randn = np.random.randn
from pandas import *

idx = Index(np.arange(1,11))
df = DataFrame(randn(10, 5), index=idx, columns=['A', 'B', 'C', 'D', 'E'])
vals = np.around(df.values,2)
#normal = colors.Normalize(vals.min()-1, vals.max()+1)
normal = colors.Normalize(0., 1.)
#normal = colors.LogNorm(vals.min()-1, vals.max()+1)

fig = plt.figure(figsize=(15,8))
ax = fig.add_subplot(111, frameon=True, xticks=[], yticks=[])

the_table=plt.table(cellText=vals, rowLabels=df.index+1, colLabels=df.columns[1:], 
                    colWidths = [0.1]*vals.shape[1], loc='center', 
                    cellColours=plt.cm.Blues(normal(vals)))

plt.show()

# df = pd.read_csv('emissions_10.txt', sep='\t')
vals = df.values[:,1:]
