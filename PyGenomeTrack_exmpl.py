

import pygenometracks.tracks as pygtk
import matplotlib.pyplot as plt

from os.path import join 

WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')

fn = join(DATADIR,'_tmp/bigwig2_X_2.5e6_3.5e6.bw')

fig, ax = plt.subplots(figsize=(16,6)) 

chrom,start,end = 'X',2700000,3100000


tk = pygtk.BigWigTrack({'file':fn})

tk.plot(ax,chrom,start,end)

plt.show()

