import numpy as np
import matplotlib.pyplot as plt

uv = np.loadtxt('polyline_uv.dat', dtype=float)

fig, ax = plt.subplots(1,2)
s = 1.05
for i in range(2):
    ax[i].plot(uv[:,2*i], uv[:,2*i+1], '.-')
    ax[i].plot(uv[0,2*i], uv[0,2*i+1], 'ro')
    ax[i].set_xlim(-s,s)
    ax[i].set_ylim(-s,s)
    ax[i].set_aspect('equal', adjustable='box')

plt.show()
