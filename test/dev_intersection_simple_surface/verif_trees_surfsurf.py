import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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

mrkp = 'x'#'.'
colp = 'r'

fig, ax = plt.subplots(nrows=1, ncols=2)



uv = np.array( [ (-1.5683766750412068E-002, -0.98431623324958784,       -1.0000000000000000,      -0.96863246649917578), 
                 (-0.69183283843488008,      -0.30816716156511997,       -1.0000000000000000,      -0.30816716156511992) ] )

for isurf in range(2):
    ax[isurf].plot( uv[:,2*isurf], uv[:,2*isurf+1], 'rx' )
    tree = np.loadtxt( 'tree_' + str(isurf+1) + '.dat' )
    for region in tree:
        draw_region2( region, ax[isurf] )
    ax[isurf].set_xlim( -1, 1 )
    ax[isurf].set_ylim( -1, 1 )    
    ax[isurf].set_aspect('equal', adjustable='box')

plt.show()

