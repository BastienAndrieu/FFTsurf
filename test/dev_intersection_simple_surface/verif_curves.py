import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


cl = ['r','b']
fig, ax = plt.subplots(nrows=1, ncols=2)

f = open('curves.dat','r')
nc = int(f.readline())
uvxyz = np.zeros((2,7))
for i in range(nc):
    for j in range(2):
        uvxyz[j] = [ float(a) for a in  f.readline().split() ]
    for isurf in range(2):
        ax[isurf].plot( uvxyz[:,2*isurf], uvxyz[:,2*isurf+1], '.-', color=cl[isurf] )
f.close()
for isurf in range(2):
    ax[isurf].set_xlim( -1, 1 )
    ax[isurf].set_ylim( -1, 1 )    
    ax[isurf].set_aspect('equal', adjustable='box')
plt.show()


