import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


uvxyz = np.loadtxt('uv_xyz.dat')

fig, ax = plt.subplots(nrows=1, ncols=2)
for isurf in range(2):
    ax[isurf].plot( uvxyz[:,2*isurf], uvxyz[:,2*isurf+1], 'rx', markersize=5 )
    ax[isurf].set_xlim( -1, 1 )
    ax[isurf].set_ylim( -1, 1 )    
    ax[isurf].set_aspect('equal', adjustable='box')
plt.show()



