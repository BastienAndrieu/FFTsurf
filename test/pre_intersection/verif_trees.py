import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches





suffixe = '_bezier' #'' #




fig, axes = plt.subplots(nrows=1, ncols=2)

for isurf in range(2):
    tree = np.loadtxt( 'tree_' + str(isurf+1) + suffixe + '.dat' )
    for region in tree:
        rect = patches.Rectangle( ( region[0], region[2] ), # lower-left corner
                                  region[1] - region[0],    # u-extent
                                  region[3] - region[2],    # v-extent
                                  fill = False,
                                  edgecolor = 'k',
                                  linewidth = 1 )
        axes[isurf].add_patch( rect )
    axes[isurf].set_aspect('equal', adjustable='box')
    axes[isurf].set_xlim( -1, 1 )
    axes[isurf].set_ylim( -1, 1 )
    axes[isurf].set_title( 'Surface #' + str(isurf+1),
                           size = 12, weight = 'bold' )
    axes[isurf].set_xticks( np.linspace(-1.0,1.0,5) )
    axes[isurf].set_yticks( np.linspace(-1.0,1.0,5) )
    
plt.show()

