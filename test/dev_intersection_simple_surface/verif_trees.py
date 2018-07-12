import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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

tuvxyz = np.loadtxt('tuv_xyz.dat')
print tuvxyz

treec = np.loadtxt('tree_c.dat')
trees = np.loadtxt('tree_s.dat')

#print treec
#print trees

mrkp = 'x'#'.'
colp = 'r'

fig, ax = plt.subplots(nrows=1, ncols=2)


if tuvxyz.size > 6:
    ax[0].plot( tuvxyz[:,0], np.zeros(tuvxyz.shape[0]), mrkp, color=colp )
elif tuvxyz.size > 0:
    ax[0].plot( tuvxyz[0], 0.0, mrkp, color=colp )
if treec.size > 2:
    for region in treec:
        draw_region1( region, ax[0] )
else:
    draw_region1( treec, ax[0] )
ax[0].set_ylim( -0.2, 0.2 )



if tuvxyz.size > 6:
    ax[1].plot( tuvxyz[:,1], tuvxyz[:,2], mrkp, color=colp )
elif tuvxyz.size > 0:
    ax[1].plot( tuvxyz[1], tuvxyz[2], mrkp, color=colp )
if trees.size > 4:
    for region in trees:
        draw_region2( region, ax[1] )
else:
    draw_region2( trees, ax[1] )
ax[1].set_xlim( -1, 1 )
ax[1].set_ylim( -1, 1 )    
ax[1].set_aspect('equal', adjustable='box')


plt.show()

