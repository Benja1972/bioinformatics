import pandas
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import os
import yaml


om1= 2*np.pi/3600. # min
om2 = om1/12.      # hour
om3 = 2*np.pi/60.  # sec

t = np.arange(0,12*3600,0.005)

phi1 = (om1*t) % (2*np.pi)
phi2 = (om2*t) % (2*np.pi)
phi3 = (om3*t) % (2*np.pi)

del1 =(phi1-phi2) % (2*np.pi)
del2 = (phi3-phi1) % (2*np.pi)

#~ plt.figure(num=2, figsize=(55, 6), dpi=80, facecolor='w', edgecolor='k')
#~ plt.plot(t,del1)
#~ plt.plot(t,del2, 'g')
#~ plt.plot(t,np.pi*, 'r')


plt.figure(num=3, figsize=(5, 6), dpi=80, facecolor='w', edgecolor='k')
plt.plot(del1,del2, '.r', MarkerSize=0.5)
plt.plot(np.pi,np.pi/2, 'ob')
plt.xlabel('$\Delta_1=\phi_{min}-\phi_{hour}$')
plt.ylabel('$\Delta_2=\phi_{sec}-\phi_{min}$')
plt.axis('square')
plt.show()


#~ plt.figure(num=1, figsize=(40, 6), dpi=80, facecolor='w', edgecolor='k')
#~ plt.plot(t,phi1)
#~ plt.plot(t,phi2, 'g')
#~ plt.plot(t,phi3, 'r')
