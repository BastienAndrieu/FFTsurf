import numpy as np
import matplotlib.pyplot as plt


xy = np.loadtxt( 'xy.dat', dtype=float )

xymin = np.amin( xy, axis=0 )
xymax = np.amax( xy, axis=0 )
rng = xymax - xymin
for i in range(2):
    if rng[i] > 2.0e-16:
        xy[:,i] = xy[:,i] / rng[i]

fig, ax = plt.subplots(1,2)

for i in range(2):
    hull = np.loadtxt( 'hull'+str(i+1)+'.dat', dtype=int ) - 1
    hull = np.hstack( ( hull, hull[0] ) )
    ax[i].plot( xy[:,0], xy[:,1], 'k.' )
    ax[i].plot( xy[hull,0], xy[hull,1], 'ro-', markerfacecolor='none' )

    xymin = np.amin( xy, axis=0 )
    xymax = np.amax( xy, axis=0 )
    rng = 0.5 * ( xymax - xymin )
    ctr = 0.5 * ( xymax + xymin )
    
    rng = 2.0 * rng
    rng = np.maximum( rng, 0.3 * np.amax( rng ) )
    xymin = ctr - rng
    xymax = ctr + rng
    
    ax[i].set_xlim( xymin[0], xymax[0] )
    ax[i].set_ylim( xymin[1], xymax[1] )
    ax[i].set_aspect('equal', adjustable='box')

plt.show()
