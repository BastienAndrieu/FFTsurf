import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D

##########################################################
def draw_region1( region, ax ):
    ax.plot( region, [0,0], '+-', color='k' )
    return
##########################################################
def draw_region2( region, ax ):
    rect = patches.Rectangle( ( region[0], region[2] ), # lower-left corner
                              region[1] - region[0],    # u-extent
                              region[3] - region[2],    # v-extent
                              fill = False,
                              edgecolor = 'k',
                              linewidth = 0.5 )
    ax.add_patch( rect )
    return
##########################################################
def draw_region3( region, ax ):
    e = np.zeros((3,2))
    for i in range(3):
        e[i,:] = region[2*i:2*(i+1)]
        for j in range(2):
            ij = (i+1)%3
            e[ij,:] = region[2*ij + j]
            for k in range(2):
                ik = (i+2)%3
                e[ik,:] = region[2*ik + k]
                ax.plot3D( *zip(e[:,0],e[:,1]), color='k', linewidth=0.5 )
    return
##########################################################







p = np.loadtxt( 'points.dat', skiprows=1 )
tree = np.loadtxt( 'tree.dat' )

if len(p.shape) < 2:
    dim = 1
else:
    dim = p.shape[1]



colp = 'r'#'b'


fig = plt.figure()
if dim == 1:
    ax = fig.gca()
    ax.plot( p, np.zeros(len(p)), '.', color=colp )
    for region in tree:
        draw_region1( region, ax )
    ax.set_ylim( -0.2, 0.2 )

if dim == 2:
    ax = fig.gca()
    ax.plot( p[:,0], p[:,1], '.', color=colp )
    for region in tree:
        draw_region2( region, ax )
        #rect = patches.Rectangle( ( region[0], region[2] ), # lower-left corner
        #                          region[1] - region[0],    # u-extent
        #                          region[3] - region[2],    # v-extent
        #                          fill = False,
        #                          edgecolor = 'k',
        #                          linewidth = 1 )
        #ax.add_patch( rect )

elif dim == 3:
    ax = fig.gca(projection='3d')
    ax.plot( p[:,0], p[:,1], p[:,2], '.', color=colp )
    for region in tree:
        draw_region3( region, ax )

ax.set_aspect('equal', adjustable='box')
plt.show()
