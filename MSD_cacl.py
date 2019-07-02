import numpy as np
import matplotlib.pyplot as plt

def msd_straight_forward(r):
    shifts = np.arange(len(r))
    msds = np.zeros(shifts.size)    

    for i, shift in enumerate(shifts):
        diffs = r[:-shift if shift else None] - r[shift:]
        sqdist = np.square(diffs).sum(axis=1)
        msds[i] = sqdist.mean()

    return msds


N = 10000
r = np.cumsum(np.random.choice([-1., 0., 1.], size=(N, 3)), axis=0)

fig = plt.figure()
plt.plot(msd_straight_forward(r))
plt.title(r'MSD($\Delta$ t)')

plt.show()
