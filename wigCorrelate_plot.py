#!/usr/bin/env python



import os
import yaml
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('ggplot')
#plt.style.use('fivethirtyeight')
#plt.style.use('seaborn-paper')

params = {'font.size': 18,
          'lines.linewidth': 3
         }
plt.rcParams.update(params)          # Plot parameters

np = '05peak/wigCorrelation_stat.txt'
bp = '06broad_peak/wigCorrelation_stat.txt'
stat_NP =  pd.read_csv(np, sep='\t', names=['repl1', 'repl2', 'corr'])
stat_BP =  pd.read_csv(bp, sep='\t', names=['repl1', 'repl2', 'corr'])

labels_NP = list(stat_NP.index) 
for i in stat_NP.index:
	labels_NP[i] =  list(stat_NP.repl1)[i][7:-12]

labels_BP = list(stat_BP.index) 
for i in stat_BP.index:
	labels_BP[i] =  list(stat_BP.repl1)[i][13:-12]


#stat_NP.plot(kind='bar')
plt.figure()
plt.plot(stat_NP.index,stat_NP['corr'],'-o', MarkerSize=10)
plt.ylabel('correlation repl1/repl2')
plt.xticks(stat_NP.index,labels_NP, rotation=75)
plt.title(list(stat_NP.repl1)[i][7:10]+' Narrow peaks')
plt.tight_layout()
plt.savefig('im/'+ list(stat_NP.repl1)[i][7:10]+'_corr_repl{1,2}'+'_NP.eps')
#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

#stat_BP.plot(kind='scatter')

#figsize=[10,8]
plt.figure()
plt.plot(stat_BP.index,stat_BP['corr'],'-o', MarkerSize=10)
plt.ylabel('correlation repl1/repl2')
plt.xticks(stat_BP.index,labels_BP, rotation=75)
plt.title(list(stat_NP.repl1)[i][7:10]+' Broads peaks')
plt.tight_layout()
plt.savefig('im/'+ list(stat_NP.repl1)[i][7:10]+'_corr_repl{1,2}'+'_BP.eps')
#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)



plt.show()




